# import imp
# import sys
# fn_, path, desc = imp.find_module('mymodule', ['/data/module/'])
# print fn_,path,desc
# mod = imp.load_module(''mymodule'', fn_, path, desc)
# print dir(mod)

#这样就会把/data/module/mymodule.py模块导入进来，load_modul方法的第一个参数可以任意写，例如mym,

#作用相当于 import mymodule as mym

#! /usr/bin/env python
# -*- coding: utf-8 -*-
# 代码有问题别找我！虽然是我写的，并且我觉得它是没问题的，如果不是你的操作原因，或许是它自己长歪了！

"""
A pipeline uses for writing NewBorn v5 vip report for each sample which base anno by BGIANNO !!!
"""

from __future__ import division
import sys, os, logging, imp, zipfile, json, re, sqlite3, shutil, errno
from optparse import OptionParser, OptionGroup
import xlsxwriter, time, fnmatch, gzip, copy, xlrd
from xml.sax.handler import ContentHandler
from xml.sax import make_parser
from collections import defaultdict
import pysam, tempfile
import coloredlogs
import haps_check,orcl_conn
coloredlogs.install()
#import multiprocessing

if sys.getdefaultencoding() != 'utf-8':
	reload(sys)
sys.setdefaultencoding('utf-8')
logger = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s: %(levelname)s: %(message)s')
logging.root.setLevel(level=logging.INFO)
logger.info("Running: python %s" % ' '.join(sys.argv))


#def cpu_count():
#	return multiprocessing.cpu_count()
drug_and_dis=['OTC','MT-RNR1','ARG1','ASS1','NAGS','CPS1','ASL','G6PD']   #药物和疾病共有基因8个
drug_not_dis=['VKORC1','HLA-B','POLG','CYP2D6','NAT2','CYP4F2','TPMT','CYP2C9','CYP2C19'] #药物非疾病基因9个

def get_alt(strings): #获取alt变异
	if re.compile('alt:(.*)\s\)').search(strings):
		m = re.compile('alt:(.*)\s\)').search(strings)
		return m.group(1)
	else:
		return strings
def extract(file_name, dirname):
	try:
		if not os.path.exists(dirname): os.makedirs(dirname)
		filein = zipfile.ZipFile(os.path.realpath(file_name))
		for f in filein.namelist():
			message = filein.extract(f, dirname)##拷贝文件，返回那个文件

		filein.close()
		return True
	except:
		return False


def parsing_conf(config_file): ##加载json
	try:
		with open(config_file) as config_files:
			config = json.load(config_files, encoding="utf-8")
		return config
	except:
		return False


def parsing_demo(module):
	try:
		mod_path, mod_name = os.path.split(module)
		mod_name, mod_ext = os.path.splitext(mod_name)
		mod_file, mod_filename, mod_data = imp.find_module(mod_name, [mod_path])
		if not mod_name in sys.modules: ##sys.modules 记录加载进来的模块，imp用来加载模块
			base = imp.load_module(mod_name, mod_file, mod_filename, mod_data)
		else:
			logger.info("[ERRO]: %s module have existed in your python, please rename your base_config !!" % mod_name)
			sys.exit(1)
		return base
	except:
		return False


def load_demo_config(demo, config, base_xml, outdir):
	demo_return = extract(demo, outdir)##拷贝demo
	config_return = parsing_conf(config)##解析
	base_xml_return = parsing_demo(base_xml)##加载必要的模块
	assert (demo_return and config_return and base_xml_return), "Erro" ##断言一定为真，否则报错
	return demo_return, config_return, base_xml_return


class NormXmlFormat(object): ##特殊字符串转化
	def __init__(self, types=True):
		adict = {
		r'&': r'&amp;',
		r'>': r'&gt;',
		r'<': r'&lt;',
		r'"': r'&quot;',
		r"'": r"&apos;"
		}
		if types is False:
			adict = dict([(v, k) for k, v in adict.iteritems()])

		self.adict = adict
		self.rx = self.make_rx()

	def make_rx(self):
		return re.compile('|'.join(map(re.escape, self.adict))) ##如果字符串很长且包含很多特殊字符，而你又不想输入一大堆反斜线，可以使用这个函数re.escape

	def one_xlat(self, match): ##match 是必要的参数吧
		return self.adict[match.group(0)]

	def __call__(self, text):
		return self.rx.sub(self.one_xlat, re.sub("(\s+)?\t(\s+)?", "\t", text)) ##为啥self.one_xlat没有参数


class ParsingDatabase(object):
	def __init__(self, data, model, token_pathogenic=set()):
		self.dbmodel = model
		self.dbdir = data
		self.mutdb = defaultdict(list)
		self.genes = defaultdict(set)
		self.diseasedb = defaultdict(list)
		self.gene_diseasedb = defaultdict(list)
		self.drug = defaultdict(dict)
		self.entertainment = defaultdict(list)
		self.flank=defaultdict(dict)
		self.token_pathogenic = token_pathogenic
		self.dbrawtitle = defaultdict(lambda: 'unkown')
		self.database_title = defaultdict(int)
		self.definition_config=defaultdict(list)
		##self.dbdir===/ifs4/HST_5B/USERS/zhaocece/Production/newborn_datasets
		mut_db = os.path.join(self.dbdir, 'nm_cds_db.txt') ##过去突变位点的详细注释
		flank_db = os.path.join(self.dbdir,'flank_db.txt')##位点侧翼序列
		disease_db = os.path.join(self.dbdir, 'disease_all_' + self.dbmodel)  ##疾病的中英文，描述NUM	name	中文名称	中文描述	参考文献	disease information	reference
		drug_db = os.path.join(self.dbdir, 'drug')
		##drug_rang是RS号位点，区间，扩展区域
		##drug_gene是基因列表，纯基因名字，有不同版本
		disease_gene = os.path.join(self.dbdir, 'disease_gene_' + self.dbmodel) #基因名列表
		drug_gene = os.path.join(self.dbdir, 'drug_gene_' + self.dbmodel)
		#drug_gene_v1是17个基因，disease_gene_v1是88个基因
		#drug_gene_v5是17个基因，disease_gene_v5是167个基因，安馨可本身药物基因是17个，与疾病166基因， overlap 8个
		gene_disease = os.path.join(self.dbdir, 'gene_disease')
		#BCKDHA	枫糖尿病	AR，基因，疾病以及遗传模式

		db_definition = os.path.join(self.dbdir, 'db_definition')  #Pathogenic	missense	$gene基因的$muntation$突变是$disease$的$type$突变。该突变为错义突变，已有该位点致病性的相关文献报道[$num$]。	$gene基因的$muntation$突变是$disease$的$type$突变。该突变为错义突变，暂未发现该位点致病性的相关文献报道。 #不同的突变类型，有文献与无文献的陈述形式
		entertainment = os.path.join(self.dbdir, 'entertainment/entertainment_last.db_' + self.dbmodel)
		self.connet = sqlite3.connect(entertainment) #连接数据库或者创建数据库
		self.con = self.connet.cursor()
		#该例程创建一个 cursor，将在 Python 数据库编程中用到。
		#数据库怎么用的还不清楚，数据库里的内容导出看看是什么，以及整合的时候用了什么字段
		norm_xml_format = NormXmlFormat(True)
		try:
			for lines in file(db_definition):
				line=lines.strip().split('\t')
				self.definition_config[(line[0].lower(),line[1].lower())] = line ##创建列表字典，根据第一列致病性字符串以及第二列突变类型
		except Exception as err:
			logger.info("[ERRO]: definition_config database format erro :%s !!" % err)

		try:
			with open (gene_disease) as f:
				for lines in f:
					line=lines.strip().split('\t')
					self.gene_diseasedb[line[0]]=line
					##BCKDHA	枫糖尿病	AR
				## 基因，疾病以及遗传模式
		except Exception as err:
			logger.info("[ERRO]: gene_disease database format erro :%s !!" % err)	

		try: 
			for lines in file(flank_db):
				line=lines.strip().split('\t')
				self.flank[line[0].lower()][line[1]]=line[2] ##染色体 位置 侧翼序列
		except Exception as err:
			logger.info("[ERRO]: flank database format erro :%s !!" % err)

		try:
			for lines in file(disease_db):
				line = norm_xml_format(lines).strip().split('\t') ##疾病中英文描述及文献，格式化输出？格式化不知道怎么实现的
				self.diseasedb[line[1].lower()] = line
		except Exception as err:
			logger.info("[ERRO]: disease database format erro :%s !!" % err)
			raise

		try:
			with open(disease_gene) as d:
				for g in d:
					g = g.strip()
					g = re.sub('^\s+|\s+$', '', g)
					self.genes['disease'].add(g)
					self.genes['all'].add(g)
			with open(drug_gene) as d:
				for g in d:
					g = g.strip()
					g = re.sub('^\s+|\s+$', '', g)
					self.genes['drug'].add(g)
					self.genes['all'].add(g)
		except  Exception as err:
			logger.info("[ERRO]: genes database format erro :%s !!" % err)
			raise

		try:
			haplotye = defaultdict(dict) ##yaowu haplotype[gene+hom][*1] or haplotype[chr start end ref call][*1]=sore转换的值
			haplotye_score = defaultdict(dict)##yaowu haplotype_score[gene][*1]
			default_gene_hap = defaultdict(lambda: ".") ##bushi yaowu de
			drug_tag = defaultdict(dict) ##yaowu drug_tag[gene][*1]=B
			genotype = defaultdict(dict)##yaowu genotype[gene][AB]=[IM,kuaidaixiexing]
			interpretation = defaultdict(dict)##yaowu interpretation[gene][BM]=[gene,ami,IM,jianshao,]
			backgroud = defaultdict(list) ##yaowu backgroud[kamaxiping]=[]
			#读取数据库，存储起来
			for line in file(os.path.join(drug_db, 'db_haplotype')): #/ifs4/HST_5B/USERS/zhaocece/Production/newborn_datasets/drug/db_haplotype
				##Gene	Chr	Start	End	Ref	Mut	NM	Haplotype	Cds	Code	Score ##Haplotype是啥意思 *2怎么来的
# CYP2C19	.	.	.	.	Ref	.	*1	.	A	0
# CYP2C19	chr10	96541615	96541616	G	A	NM_000769.1	*2	c.681G>A	B	1
# CYP2C19	chr10	96540409	96540410	G	A	NM_000769.1	*3	c.636G>A	B	1

				if line.startswith("#"):
					continue
				lines = re.sub("(\s+)?\t(\s+)?", "\t", line).strip().split('\t')
				gene, chrom, start, end, ref, mut, trans, haplo, chgvs, tag, score = lines
				if mut == "Ref":
					default_gene_hap[gene] = haplo
				haplotye_score[gene][haplo] = float(score) ##CYP2C19  *3 =1
				drug_tag[gene][haplo] = tag ##CYP2C19  *3 =B
				if chrom is not "." and start is not "." and end is not ".":
					tmp_mut = ":".join([chrom, start, end, ref, mut])
					haplotye[tmp_mut][haplo] = 1.0 / haplotye_score[gene][haplo] if haplotye_score[gene][haplo] > 0.0 \
						else haplotye_score[gene][haplo]
				elif haplo == "Hom":
					tmp_mut = ":".join([gene, "纯合"])
					haplotye[tmp_mut][haplo] = 1.0 / haplotye_score[gene][haplo] if haplotye_score[gene][haplo] > 0.0 \
						else haplotye_score[gene][haplo]

			for line in file(os.path.join(drug_db, 'db_genotype')):
				# Gene	Code	Genotype	Description
				# CYP2C19				AA				EM				快代谢型
				# CYP2C19				AB				IM				中间代谢型
				# CYP2C19				BA				IM				中间代谢型
				# CYP2C19				AD				Unknown				临床意义未明
				# CYP2C19				DA				Unknown				临床意义未明
				# CYP2C19				BB				PM				弱代谢型
				if line.startswith("#"):
					continue
				lines = re.sub("(\s+)?\t(\s+)?", "\t", line).strip().split('\t')
				gene, code, genotypes, description = lines
				genotype[gene][code] = [genotypes, description]

			for line in file(os.path.join(drug_db, 'db_interpretation')):
				##Gene	Drug	Genotype	Short recommendation ch	Therapeutic (dose) recommendation	Reference	Short recommendation	Therapeutic (dose) recommendation	Reference
#CYP2D6	阿米替林	IM	减少剂量	CPIC（Clinical Pharmacogenetics Implementation Consortium）剂量指导CYP2D6中间代谢型使用阿米替林时，建议降低25%的起始剂量，并使用治疗药物指导剂量调整。与CYP2D6抑制剂药物联合使用时，需要适当降低阿米替林的剂量[1,2]。	[1]Hicks JK, Swen JJ, Thorn CF et al. Clinical Pharmacogenetics Implementation Consortium guideline for CYP2D6 and CYP2C19 genotypes and dosing of tricyclic antidepressants. Clin Pharmacol Ther 2013; 93: 402-408.[2]Par Pharmaceutical, Inc.  Chlordiazepoxide and Amitriptyline drug label.2011.	Reduce dose	The CPIC Dosing Guideline for amitriptyline a 25% dose reduction and monitor plasma concentration or select alternative drug (e.g., citalopram, sertraline) for CYP2D6 intermediate metabolizers[1,2].	[1]Hicks JK, Swen JJ, Thorn CF et al. Clinical Pharmacogenetics Implementation Consortium guideline for CYP2D6 and CYP2C19 genotypes and dosing of tricyclic antidepressants. Clin Pharmacol Ther 2013; 93: 402-408. [2]Par Pharmaceutical, Inc.  Chlordiazepoxide and Amitriptyline drug label.2011.


				if line.startswith("#"):
					continue
				rows = norm_xml_format(line).strip().split('\t')
				gene, drug_ch, drug_gen, short_ch, long_ch, reference_ch, short, long, reference = rows

				if drug_gen not in interpretation[gene]:
					interpretation[gene][drug_gen] = [rows]
					#将基因和基因型（不是变异的那个基因型，而是跟用药有关的基因型）为键，值为所有的解读信息的list
				else:
					interpretation[gene][drug_gen].append(rows) ##如果已经存在则加上这个list

			for line in file(os.path.join(drug_db, 'db_background')):
				#药物信息
				##序号	药物	Drug	背景	参考文献	Drug info	Reference
#1	卡马西平	Carbamazepine	卡马西平属于抗惊厥和抗癫痫药物，适用于部分性和全身性癫痫发作，同时也可作为治疗缓解三叉神经痛和舌咽神经痛的特效药物，另外用于预防或治疗躁狂-抑郁症、中枢性部分尿崩症等症状。由卡马西平引起的常见不良反应包括头晕、共济失调、嗜睡和疲劳等症状。另外有极少部分病人会产生严重、甚至致死的变态反应，包括史蒂文斯-约翰逊综合征（SJS）和中毒性表皮坏死溶解症(TEN)。2004年首次在汉族中发现由卡马西平引起的SJS/TEN与HLA-B*1502基因型紧密关联[1]，其后同一个研究小组后续的研究以及在香港汉族人中研究进一步验证了这个结论[2,3]。2011年在台湾招募4877名未服用卡马西平的候选受试者，经过HLA*1502检测分型后，建议阳性携带者不服用卡马西平而是服用备选的另外一种药物，而阴性受试者继续服用卡马西平，结果显示没有出现一例SJS/TEN患者，而由卡马西平导致的SJS/TEN历史发病率为0.23%[4]。这项研究表明了对于汉族人使用卡马西平前进行HLA*1502基因型分型的必要性。

				if line.startswith("#"):
					continue
				rows = norm_xml_format(line).strip().split('\t')
				backgroud[rows[1]] = rows

			self.drug['haplotye'] = haplotye
			self.drug['haplotye_score'] = haplotye_score
			self.drug['drug_tag'] = drug_tag
			self.drug['genotype'] = genotype
			self.drug['interpretation'] = interpretation
			self.drug['backgroud'] = backgroud
			self.drug['default_gene_hap'] = default_gene_hap
		except Exception as err:
			logger.info("[ERRO]: Drug database format erro :%s !!" % err)
			raise

		try:
			enter_sql = "SELECT * FROM commits"
			self.con.execute(enter_sql)
			enter_meassage = self.con.fetchall()
			for rows in enter_meassage:
				#('TERMID_1', 'chr17', '42328620', '42328621', 'G', 'A', 'Hom', 'SLC4A1', 'rs2285644', 'AA', '2', '0.0')
				num, chrom, start, end, ref, mut, zygosity, gene, rs, g_type, result, fre = rows
				key = (chrom, start, end, ref, mut, zygosity)
				self.entertainment[key].append(rows) ##这个数据库不知道从哪来的，个体特征？存储突变的位置refalt突变以及纯杂合，值为整行最后两列不太确定，尤其倒第二列
		except Exception as err:
			logger.info("[ERRO]: entertainment database select erro :%s !!" % err)
			raise

		try:
			rawtitle = defaultdict(str)
			database_title = defaultdict(int)
			tmp_mut_db = [line for line in file(mut_db)] ##过去突变位点的详细注释
			mut_db_title = norm_xml_format(tmp_mut_db.pop(0)).strip().split('\t') #biaotou格式化
			for number_of_terms in xrange(len(mut_db_title)):
				t_name = mut_db_title[number_of_terms].lower()
				rawtitle[t_name] = mut_db_title[number_of_terms]
				database_title[t_name] = number_of_terms
			self.database_title = database_title ##过去突变位点的详细注释的表头的对应哪个列的记录
			self.dbrawtitle = rawtitle  ##现在的小写后的表头对应的过去突变位点的详细注释的表头
			for lines in tmp_mut_db: ##一下是每一列的格式规整，提取对应字段
				line = norm_xml_format(lines).strip().split('\t')
				chrom = 'chr' + re.sub('^chr', '', line[database_title['chrom']])
				start = line[database_title['start']]
				end = line[database_title['stop']]
				trans = re.sub('\.[0-9]+$', '', line[database_title['nbtrans']])
				phgvs = re.sub('\.[0-9]+$', '', line[database_title['transcript']])
				ref = line[database_title['ref']]
				call = line[database_title['call']]
				definition = line[database_title['definition']]  ##致病性
				gene_sym = line[database_title['#gene']]
				if len(token_pathogenic) and definition.lower() not in self.token_pathogenic:  ##不知道self.token_pathogenic怎么获取到的，致病性的种类
					continue
				if gene_sym not in self.genes['all']: #过滤掉非目的基因的
					continue
				term_keys = (chrom, start, end, ref, call)
				if (term_keys not in self.mutdb) or (trans == phgvs): ##如果条目不在已有的数据库或者两列的转录本一致，存储成数据库mutdb，键chrom, start, end, ref, call
					self.mutdb[term_keys] = line
		except Exception as err:
			logger.info("[ERRO]: Mut database format erro :%s !!" % err)
			raise

	def __del__(self):
		self.connet.close() ##关闭数据库的关联

	def description(self, num, factor):
		sql = "SELECT * FROM descriptions WHERE Num = \"" + num + "\" AND Factor = \"" + factor + "\""
		try:
			self.con.execute(sql)
			results = self.con.fetchone()
			return "\t".join(results)
		except Exception as err:
			logger.info("[ERRO]: Commits select erro : %s" % err)
			return None

def gender_message(sample): #获取性别
	gender_dic={}
	for rows in file('final_result'):
		#Order   Sample  Ave_depth       rmdup_depth     Cover   Cover_30X       Gender  Index   Mass_result
# 1       18S0150055R2    852.39  202.07  99.71%  98.41%  F       18L0150055R2-113        NA
# 2       18S0415555      761.26  194.55  99.25%  96.29%  F       18L0415555-108  NA
# 3       18S0409209      613.46  201.79  99.55%  97.31%  F       18L0409209-118  NA

		if rows.startswith('Order'):continue
		row = rows.strip().split('\t')
		gender_dic[row[1]] = row[6][-1]
	return gender_dic[sample]
	
def sample_messages(sample_list, message_file=None, title=None): #ganmade
	if not title:
		title = dict()
	norm_xml_format = NormXmlFormat(True)
	back_xml_format = NormXmlFormat(False)
	message = defaultdict(dict)
	messages = ['born_time','guardian_name','name']
	connect = orcl_conn.Oracle() ##实例
	with open(sample_list) as samples:##读取list时候，抽提客户的记录信息
		for sample in samples:
			rows = sample.strip().split('\t')
			sample_name = rows[0]
			sql = "select * FROM DX_SIMS.VIEW_EAR_REPORT_ANALYSIS where MAIN_SAMPLE_NUM = \'"+sample_name+"\'"
			vip_info = connect.vip_result(sql) ##没有测试出来
			gender = rows[1]
			#name=rows[5].decode('gb2312')
			message[sample_name] = defaultdict(lambda: "-")
			message[sample_name]['gender'] = "女" if re.compile('f', re.I).match(gender) else '男'
#			message[sample_name]['guardian_name'] = "缺省" if re.compile('null', re.I).match(guardian_name) else guardian_name
			for i in range (len(messages)):
				message[sample_name][messages[i]] = vip_info[i] ##orcl数据库里有客人的出生日期，名字以及guardian_name
			#message[sample_name]['name'] = "-" if re.compile('null', re.I).match(name) else name
	if message_file: ##这是个什么文件
		with open(message_file) as messages:
			for line in messages:
				rows = norm_xml_format(line).strip().split('\t')
				sample = rows[title['sample_name']] if 'sample_name' in title else rows[0]
				sample = back_xml_format(sample)
				for k, v in title.iteritems():
					message[sample][k] = rows[v]
	return message


def xlsx2list(filein, sheet_name=set()): #读取sheet_name指定的sheet，展开成列表
	def check_value_from_xlrd(values):
		if type(values) is float and values == int(values):
			values = int(values)
		return str(values)

	tables = defaultdict(dict)
	if os.path.isfile(filein + '.xls'):
		filein = filein + '.xls'
	elif os.path.isfile(filein + '.xlsx'):
		filein = filein + '.xlsx'
	if not filein.endswith('xlsx') and not filein.endswith('xls'):
		return tables
	data = xlrd.open_workbook(filein)
	for sheet in data.sheets():
		if len(sheet_name) and sheet.name not in sheet_name:
			continue
		nrows = sheet.nrows
		ncols = sheet.ncols
		if nrows < 1 or ncols < 1:
			continue
		tables[sheet.name] = [[check_value_from_xlrd(sheet.cell(i, j).value).replace('\r', '').replace('\n', '')
		                       for j in range(ncols)] for i in range(nrows)] ##行列循环，修正值，并且展开成一个列表
	return tables


def sample_gene_id(geneid=None): ##读取id.xls的内容，记录成每个样本在每个基因（或者基因组合）的基因型？比如g-T-A g-C-G，并返回
	gene_id_dict = defaultdict(dict)
	if not geneid:
		return gene_id_dict
	norm_xml_format = NormXmlFormat(True)
	sample_gene_id_messages = xlsx2list(geneid, {"haps"})##xlsx2list读取sheet_name指定的sheet，展开成列表
	sample_gene_id_message = sample_gene_id_messages["haps"]
	if len(sample_gene_id_message):
		gene_id_title = sample_gene_id_message.pop(0)
		if gene_id_title.pop(0) == "title":
			for rows in sample_gene_id_message:
				sample_name = rows.pop(0)
				if len(rows) != len(gene_id_title):
					continue
				for num in xrange(len(rows)):
					word = '-'.join(gene_id_title[num].split()).upper()
					gene_id_dict[sample_name][word] = norm_xml_format(rows[num].upper()).strip()
		else:
			logger.info("[ERRO]: Gene ID file format erro !!")
	else:
		logger.info("[ERRO]: Gene ID file format erro !!")
	return gene_id_dict


def findfile(pathdir, sample, name, type, filters=None):
	for root, dirs, files in os.walk(os.path.join(pathdir, sample)):
		for filename in fnmatch.filter(files, str(name) + '*.' + type):
			if os.path.isdir(filename):
				continue
			if filters and re.compile(filters + '.' + type + "$", re.I).search(filename):
				continue
			return os.path.join(root, filename)
	return None


def indexfile(variation, bamdir, sample_name, indexdir):
	vcffile, depthfile = findfile(variation, sample_name, str(sample_name) + '.', 'gz', 'vcf'), \
	                     findfile(bamdir, sample_name, str(sample_name) + '.', 'bam')
	if vcffile is None or depthfile is None:
		return vcffile, depthfile
	if not os.path.isfile(depthfile + '.bai') and not os.path.isfile(re.sub("bam$", 'bai', depthfile)):
		filedir, input_file = os.path.split(depthfile)
		if not os.access(filedir, os.W_OK):#os.W_OK 包含在access()的mode参数中 ， 测试path是否可写。
			os.symlink(depthfile, os.path.join(indexdir, input_file))
			depthfile = os.path.join(indexdir, input_file)
		retval = pysam.index(depthfile)##不确定这个depthfile是啥，就是bam文件
	return vcffile, depthfile


def copy_file(olddir, new_dir):
	tmp_time = time.strftime('%Y.%m.%d.%H.%M', time.localtime(time.time()))
	try:
		if os.path.exists(new_dir):
			if os.path.exists("_".join([new_dir, 'bak', tmp_time])):
				shutil.rmtree("_".join([new_dir, 'bak', tmp_time]), ignore_errors=True)
			os.rename(new_dir, "_".join([new_dir, 'bak', tmp_time]))
		shutil.copytree(olddir, new_dir)
	except OSError as err:
		if err.errno == errno.ENOTDIR:
			shutil.copy(olddir, new_dir)
		else:
			raise


def isnum(v):
	try:
		x = float(v)
	except TypeError:
		return False
	except ValueError:
		return False
	except Exception, e:
		return False
	else:
		return True

def parsing_mut_indb(batch_num, sample, indb, check_gene, gender, newborndb, haplotype,  , f_out, all_indb_out):
	##batch_num 批次或者项目编号
	##sample 样品名称
	##indb是个临时汇总记录的字典，如果不在外部数据库或者转录本是目的转录本，则记录
	####check_gene在已知库中的位点，然后基因与疾病组合键的包含的位点list集合
	##message['gender'] 性别信息
	##newborndb是流程数据库路径下的各种信息的组合
	##haplotype 项目分析的result路径下有id.xls 以及每个样本的好像是基因型但是就两行
	##titles bedanno的表头，直接看起来就是，应该没有其他改动
	##tmp_out
	##all_indb_out在已有的数据库中的记录
	fout = open(f_out, 'w')
	annotations = defaultdict(set)
	drug_gene_hap = defaultdict(lambda: ".")
	drug_haplo_score_a = defaultdict(float)
	drug_haplo_score_b = defaultdict(float)##存？
	other_drug_gene = set()
	drug_usage = defaultdict(set)
	default_hap = newborndb.drug['default_gene_hap'] ##如果没突变使用默认的状态
	check_nat2 = True
	back_xml_format = NormXmlFormat(False)
	pathogenic_trans = {
	"pathogenic": "已知致病",
	"dm": "已知致病",
	"likely pathogenic": "疑似致病",
	"dp": "疑似致病",
	"drug": "药物"
	}
	database_title = sorted(newborndb.database_title.keys(), key=lambda s: int(newborndb.database_title[s]))##过去的记录，已经解读过的位点的详细信息的表头，已有数据库表头是标杆
	dis_genes = newborndb.genes['disease']##166gene
	drug_genes = newborndb.genes['drug'] ##drug基因，drug_gene_v5 17个
#	fout.write(back_xml_format("\t".join([newborndb.dbrawtitle[i] for i in database_title])) + '\n')
	fout.write("Project" + '\t' + "Sample" + '\t' + back_xml_format("\t".join([newborndb.dbrawtitle[i] for i in database_title if newborndb.database_title[i] <= newborndb.database_title["phgvs"]]))) ##转换成原来的格式
	fout.write('\t'+"ForPrimerDesign" + '\t' + back_xml_format("\t".join([newborndb.dbrawtitle[i] for i in database_title if newborndb.database_title[i] >newborndb.database_title["phgvs"]])) + '\n')
	for annotation_key, annotation_value in indb.iteritems():##indb是个临时汇总记录的字典，如果突变位点在外部数据库或者再加上转录本是目的转录本，则记录,键是chr start end ref alt
		database_line = newborndb.mutdb[annotation_key]##仍然是和indb感觉一样呢，##如果条目不在已有的数据库或者两列的转录本一致，存储成数据库mutdb，键chrom, start, end, ref, call，indb是目前这个样本在数据库中的位点，mutdb是过去数据库中已经解读的位点
		#循环变异条目，假如目前解读的位点不在过去的数据库里怎么处理？
		for db_item_value, db_item_number in newborndb.database_title.iteritems(): ##过去的记录，已经解读过的位点的详细信息的表头，用过去的表头循环现在的信息，上面最外层是现在样本的突变信息
			if db_item_value in annotation_titles and not db_item_value.lower() == "flank":
				database_line[db_item_number] = annotation_value[annotation_titles[db_item_value]] ##annotation_value是个list，我目前觉得多此一举的感觉，可能某个转换没意识到
		database_line[newborndb.database_title["mutfunction"]] = database_line[newborndb.database_title["nbtrans"]]
		database_line[newborndb.database_title["nbchgvs"]] = annotation_value[annotation_titles['chgvs']]
		database_line[newborndb.database_title["nbtrans"]] = annotation_value[annotation_titles['transcript']]
		chgvs = database_line[newborndb.database_title["nbchgvs"]]
		trans = database_line[newborndb.database_title["nbtrans"]]
		term_phgvs = "."
		try:
			phgvs = database_line[newborndb.database_title["phgvs"]]
			for pv in re.findall('p.[^\s+|:]+', phgvs, re.I):  #p.Y240C | p.Tyr240Cys
				if len(pv) >= len(term_phgvs): term_phgvs = pv  ##最终取非缩写的p点？看着是，看嘛去
			phgvs = term_phgvs
		except:
			phgvs = term_phgvs
		gene_name = database_line[newborndb.database_title["#gene"]]
		inheritance = database_line[newborndb.database_title["tag"]]
		if gene_name in newborndb.gene_diseasedb : inheritance = newborndb.gene_diseasedb[gene_name][2]##遗传模式，AR AD
		gtype = database_line[newborndb.database_title["zygosity"]]
		zygosity = database_line[newborndb.database_title["zygosity"]]
		functions = database_line[newborndb.database_title["functionname"]]
		funcregion = database_line[newborndb.database_title["funcregion"]]
		flank = database_line[newborndb.database_title["flank"]] ##过去数据库中的情况
		if re.compile('hom', re.I).match(gtype):
			gtype = "纯合"
		elif re.compile('het', re.I).match(gtype):
			gtype = "杂合"
#		else:
#			continue
		if database_line[newborndb.database_title["definition"]].lower() not in pathogenic_trans:
			continue ##pathogenic dm likely pathogenic dp drug
		pathogenic = pathogenic_trans[database_line[newborndb.database_title["definition"]].lower()]
		mut_definiton = database_line[newborndb.database_title["chlongdescription"]]
		references = database_line[newborndb.database_title["reference"]]
		forPrimerDesign = (gene_name,trans,chgvs,phgvs,funcregion,zygosity,flank)
		if gene_name in dis_genes:
			disease_en = database_line[newborndb.database_title["disease"]]
			if disease_en.lower() in newborndb.diseasedb:
				disease_raw = newborndb.diseasedb[disease_en.lower()]
			else:
				disease_raw = ["."] * 7
			disease_ch = disease_raw[2] if disease_raw[2] is not "." else disease_en
			disease_descri_ch = disease_raw[3]
			disease_refer_ch = disease_raw[4]
			fout.write(batch_num + '\t' + sample + '\t')
			all_indb_out.write(batch_num + '\t' + sample + '\t')
			for i in range(len(database_line)):
				if(i == newborndb.database_title["phgvs"] ):
					fout.write(back_xml_format(database_line[i]) + '\t' + back_xml_format("; ".join(forPrimerDesign)) + '\t')
					all_indb_out.write(back_xml_format(database_line[i]) + '\t' + back_xml_format("; ".join(forPrimerDesign)) + '\t')
				else:
					fout.write(back_xml_format(database_line[i]) + '\t')
					all_indb_out.write(back_xml_format(database_line[i]) + '\t')
			fout.write('\n')
			all_indb_out.write('\n')
####		fout.write(back_xml_format("\t".join(database_line)) + '\n')
			annotations["mut_definiton"].add((disease_ch, disease_en, mut_definiton, references))
			annotations["disease_message"].add((disease_ch, disease_en, disease_descri_ch, disease_refer_ch))
			check_gene_num = 0
			for ann_key in check_gene[(gene_name, disease_en.lower())]:
				if newborndb.mutdb[ann_key][newborndb.database_title["definition"]].lower() in pathogenic_trans:
					check_gene_num += 1
			if ( check_gene_num > 1) or \
					(not re.compile('[ax]r', re.I).search(inheritance)) or \
					(re.compile('x', re.I).search(inheritance) and gender == "男") or \
					(gtype == "纯合" or gtype == "Hemi"):
				annotations["pathgenic"].add(
					(disease_ch, disease_en, gene_name, inheritance.upper(), chgvs, phgvs, functions, gtype, pathogenic,trans,references))
			else:
				annotations["carrier"].add(
					(disease_ch, disease_en, gene_name, inheritance.upper(), chgvs, phgvs, functions, gtype, pathogenic,trans,references))
		if gene_name in drug_genes:
			if (not gene_name in newborndb.drug['haplotye']) and (
			not ":".join(annotation_key) in newborndb.drug['haplotye']):
				continue
			fout.write(batch_num + '\t' + sample + '\t')
			all_indb_out.write(batch_num + '\t' + sample + '\t')
			for i in range(len(database_line)):
				if(i == newborndb.database_title["phgvs"]):
					fout.write(back_xml_format(database_line[i]) + '\t' + back_xml_format("; ".join(forPrimerDesign)) + '\t')
					all_indb_out.write(back_xml_format(database_line[i]) + '\t' + back_xml_format("; ".join(forPrimerDesign)) + '\t')
				else:
					fout.write(back_xml_format(database_line[i]) + '\t')
					all_indb_out.write(back_xml_format(database_line[i]) + '\t')
			fout.write('\n')
			all_indb_out.write('\n')
#			fout.write(back_xml_format("\t".join(database_line)) + '\n')
			mut_key = ":".join(annotation_key)
			if gene_name == "NAT2":
				check_nat2 = False
			if mut_key in newborndb.drug['haplotye']:
				for haplo, score in newborndb.drug['haplotye'][mut_key].iteritems():
					drug_haplo_score_a[(gene_name, haplo)] += float(score)
					if gtype == "纯合" or gtype == "Hemi":
						drug_haplo_score_b[(gene_name, haplo)] += float(score)
			elif ":".join([gene_name, gtype]) in newborndb.drug['haplotye']:
				for haplo, score in newborndb.drug['haplotye'][":".join([gene_name, gtype])].iteritems():
					drug_haplo_score_a[(gene_name, haplo)] += float(score)
					drug_haplo_score_b[(gene_name, haplo)] += float(score)
			else:
				other_drug_gene.add(gene_name)
	fout.close()
	drug_lane1 = defaultdict(set)
	drug_lane2 = defaultdict(set)

	for gene, hap in [key for key, value in drug_haplo_score_a.iteritems() if round(value, 3) == 1.0]:
		drug_lane1[gene].add(hap)
		if round(drug_haplo_score_b[(gene, hap)], 3) == 1.0:
			drug_lane2[gene].add(hap)
	for gene in drug_lane1:
		haps = sorted(drug_lane1[gene], key=lambda s: newborndb.drug['haplotye_score'][gene][s], reverse=True)
		if len(haps) > 1:
			drug_gene_hap[gene] = haps[:2]
		elif len(haps) and len(drug_lane2[gene]):
			drug_gene_hap[gene] = [haps[0], drug_lane2[gene].pop()]
		elif len(haps):
			drug_gene_hap[gene] = [default_hap[gene], haps[0]]
		else:
			drug_gene_hap[gene] = [default_hap[gene]] * 2
	for gene in other_drug_gene:
		if gene not in drug_gene_hap :
			drug_gene_hap[gene] = [default_hap[gene]] * 2
	if haplotype:
		sample_haplotype_messages = xlsx2list(haplotype, {"hla-abo"})
		sample_haplotype_message = sample_haplotype_messages["hla-abo"]
		if len(sample_haplotype_message):
			for rows in sample_haplotype_message:
				type1, type2, score = rows[:3]
				if not isnum(score):
					continue
				t_1 = ":".join(type1.split(":")[:2]).upper()
				t_2 = ":".join(type2.split(":")[:2]).upper()
				if t_1 in newborndb.drug['interpretation']['HLA-B']:
					test_hla = t_1
				elif t_2 in newborndb.drug['interpretation']['HLA-B'] and float(score) < 80:
					test_hla = t_2
				else:
					test_hla = t_1 if float(score) >= 50 else t_2
				if not len(test_hla):
					test_hla = "NOVEL"
				if test_hla in newborndb.drug['interpretation']['HLA-B']:
					for rows in newborndb.drug['interpretation']['HLA-B'][test_hla]:
						drug_gene, drug_ch, drug_gen, drug_short_ch, drug_long_ch, drug_reference_ch, \
						drug_short, drug_long, drug_reference = rows

						drug_en = newborndb.drug['backgroud'][drug_ch][2]
						drug_back_ch = newborndb.drug['backgroud'][drug_ch][3]
						drug_back_ref_ch = newborndb.drug['backgroud'][drug_ch][4]
						annotations["drug"].add((drug_ch, drug_en, drug_gene,'HLA-'+test_hla, '-', drug_short_ch))
						annotations["drug_definiton"].add((drug_ch, drug_en, drug_long_ch, drug_reference_ch))
						annotations["drug_message"].add((drug_ch, drug_en, drug_back_ch, drug_back_ref_ch))
						drug_usage[drug_en].add(drug_short_ch)
	for gene, hap in drug_gene_hap.iteritems():
		tag = "".join([newborndb.drug['drug_tag'][gene][t] for t in sorted(hap)])
		if hap[0] == hap[1]:
			if hap[0].startswith("*"):
				haps = "".join(hap)
			else:
				haps = hap[0]
		else:
			if hap[0].startswith("*"):
				haps = "".join(hap)
			else:
				haps = sorted(hap)[1]
		try:
			drug_genotype, drug_description = newborndb.drug['genotype'][gene][tag]
		except:
			continue
		for rows in newborndb.drug['interpretation'][gene][drug_genotype]:
			drug_gene, drug_ch, drug_gen, drug_short_ch, drug_long_ch, drug_reference_ch, \
			drug_short, drug_long, drug_reference = rows
			if re.compile('standard dose',re.I).search(drug_short): continue
			drug_en = newborndb.drug['backgroud'][drug_ch][2]
			drug_back_ch = newborndb.drug['backgroud'][drug_ch][3]
			drug_back_ref_ch = newborndb.drug['backgroud'][drug_ch][4]
			annotations["drug"].add((drug_ch, drug_en, drug_gene, haps, drug_description, drug_short_ch))
			annotations["drug_definiton"].add((drug_ch, drug_en, drug_long_ch, drug_reference_ch))
			annotations["drug_message"].add((drug_ch, drug_en, drug_back_ch, drug_back_ref_ch))
			drug_usage[drug_en].add(drug_short_ch)
	if check_nat2 :
		drug_ch = "异烟肼"
		drug_en = "Isoniazid"
		drug_gene = "NAT2"
		haps = "*4*4"
		drug_description = "快代谢"
		drug_short_ch = "增加剂量"
		drug_long_ch = "快乙酰化代谢型携带者其异烟肼清除率明显高于慢乙酰化携带者，同时对异烟肼的敏感性降低，更容易导致治疗失败；" \
		               "而慢乙酰化在常规剂量治疗中则更容易导致肝脏或神经系统毒性[1,2]。研究显示，对于NAT2快代谢型需要增加50%" \
		               "左右的异烟肼剂量，以达到预期的治疗效果[3]。"
		drug_reference_ch = "[1]Chen B, et al. The influence of NAT2 genotypes on the plasma concentration of isoniazid " \
		                    "and acetylisoniazid in Chinese pulmonary tuberculosis patients. Clin Chim Acta. 2006 Mar; 365" \
		                    "(1-2):104-8. [2]Kim SH, et al. Genetic polymorphisms of drug-metabolizing enzymes and anti-TB " \
		                    "drug-induced hepatitis. Pharmacogenomics. 2009 Nov; 10(11):1767-79. [3]Kinzig-Schippers M, et " \
		                    "al. Should we use N-acetyltransferase type 2 genotyping to personalize isoniazid doses? Antimic" \
		                    "rob Agents Chemother. 2005 May; 49(5):1733-8."
		drug_back_ch = "异烟肼是各器官系统、各类型结核病和结核病预防治疗的首选药物，适用于初治和复治的各型肺结核、结核性胸膜炎、腹膜炎" \
		               "、心包炎、及消化道、泌尿生殖器结核，骨关节结核、淋巴结结核等。是结核性脑膜炎的必选药物。本品单用易产生耐药故需" \
		               "联合其他抗结核病药物使用。常联合利福平、吡嗪酰胺和链霉素（或乙胺丁醇)治疗初治病例，亦可联合除上述药物以外的其" \
		               "他抗结核药物用于复治病例的治疗。异烟肼是治疗结核病的一线药物，最为常见的不良反应是非常容易导致神经系统毒性肝毒" \
		               "性。研究显示导致的肝毒性的主要因素与异烟肼的主要代谢酶—乙酰转移酶Ⅱ（NAT2）的基因多态性有关。根据NAT2代谢能力" \
		               "的不同分为快乙酰化和慢乙酰化代谢型，携带一个（中间代谢型）或两个（快代谢）高活性NAT2等位基因称为快乙酰化代谢型，" \
		               "携带两个低活性NAT2等位基因的称为慢乙酰化代谢型（慢代谢）[1,2]。快乙酰化代谢型携带者其异烟肼清除率明显高于慢乙" \
		               "酰化携带者，同时对异烟肼的敏感性降低，更容易导致治疗失败；而慢乙酰化在常规剂量治疗中则更容易导致肝脏或神经系统" \
		               "毒性[3,4]。"
		drug_back_ref_ch = "[1]Sanofi-aventis.Rifampin, isoniazid and pyrazinamide drug label.2010. [2]Kinzig-S" \
		                   "chippers M, et al. Should we use N-acetyltransferase type 2 genotyping to personalize isoniazid " \
		                   "doses? Antimicrob Agents Chemother. 2005 May; 49(5):1733-8. [3]Chen B, et al. The influence of NA" \
		                   "T2 genotypes on the plasma concentration of isoniazid and acetylisoniazid in Chinese pulmonary tub" \
		                   "erculosis patients. Clin Chim Acta. 2006 Mar; 365(1-2):104-8. [4]Kim SH, et al. Genetic polymorph" \
		                   "isms of drug-metabolizing enzymes and anti-TB drug-induced hepatitis. Pharmacogenomics. 2009 Nov;" \
		                   " 10(11):1767-79."
		annotations["drug"].add((drug_ch, drug_en, drug_gene, haps, drug_description, drug_short_ch))
		annotations["drug_definiton"].add((drug_ch, drug_en, drug_long_ch, drug_reference_ch))
		annotations["drug_message"].add((drug_ch, drug_en, drug_back_ch, drug_back_ref_ch))
		drug_usage[drug_en].add(drug_short_ch)
	return annotations, drug_usage


def float_to_per(rate):
	return "%.1f%%"  %(float(rate)*100)


def parsefile(file):
	parser = make_parser()
	parser.setContentHandler(ContentHandler( ))
	parser.parse(file)


def compile_dir(dir, out):
	file_out = zipfile.ZipFile(out, 'w')
	for root, dirs, files in os.walk(dir):
		for f in files:
			abs_path = os.path.join(root, f)
			rel_path = os.path.relpath(abs_path,dir)
			file_out.write(abs_path, rel_path, zipfile.ZIP_STORED)
	shutil.rmtree(dir, ignore_errors=True)
	file_out.close()


def test2xlsx(filname, *args):
	excel_name = xlsxwriter.Workbook(str(filname) + '.xlsx')
	sheet = 0
	format=excel_name.add_format({"font":'Arial', "size":11, "bold":1, "align":'left', "valign":'top',"text_wrap":0})
	format.set_border(1)
	format.set_bg_color('#FFCC99')
	format_black=excel_name.add_format({"font":'Arial', "size":11, "bold":1, "align":'left', "valign":'top',"text_wrap":0})
	format_black.set_border(1)
	format_black.set_bg_color('#CCFFFF')
	for files, formact_dict in args:
		path, name = os.path.split(files)
		names = name.split('_')
		name = "_".join(names[1:])
		if len(name) > 20:
			sheet += 1
			name = 'sheet' + str(sheet)
		worksheet = excel_name.add_worksheet(name)
		row = 0
		for line in file(files):
			rows = line.strip().split('\t')
			row += 1
			row_name = 'A' + str(row)
			if row % 2:
				worksheet.write_row(row_name,rows,format_black)
			else:
				worksheet.write_row(row_name,rows,format)
		if formact_dict is not None:
			for colnumber, (weight, otherdict) in formact_dict.iteritems():
				worksheet.set_column(colnumber, colnumber, weight, None, otherdict)
		os.remove(files)
	excel_name.close()

def auto_report(batch_num, sample, variation, bamdir, tmpdir, outdir, newborndb, haplotype, message, gene_id, config, base, attachment, report_type, varsheet_dict, varsheet_title, all_indb_out, all_Notindb_out):
	logger.info("|_ _ _ _ Start to writing reports of %s " % sample)
	copy_file(tmpdir, os.path.join(outdir, sample)) ##不知道tmpdir有什么内容
	rep = open(os.path.join(outdir, sample, 'word/document.xml'), 'w')
	attach = copy.copy(attachment)
	try:
		norm_xml_format = NormXmlFormat(True)
		back_xml_format = NormXmlFormat(False)
		report_times = time.strftime('%Y年%m月%d日', time.localtime(time.time()))
		annofile, depthfile = indexfile(variation, bamdir, sample, tmpdir)
		raw_result_xlsx = findfile(variation, sample, str(sample) + '_', 'xlsm') ##bedanno注释后的文件吧
		if raw_result_xlsx is None:
			raw_result_xlsx = findfile(variation, sample, str(sample) + '_', 'xlsx')
		if raw_result_xlsx is not None:
			attach.add(raw_result_xlsx)
		try:
			depth_from_bam = pysam.AlignmentFile(depthfile, 'rb')
		except Exception, e:
			logger.info("|_ _ _ _ Bam file error : %s - %s " % (sample, e))
			depth_from_bam = None
		annotations = defaultdict(list)
		partlast = config["disease_drug_table"]
		last_table = defaultdict(bool)
		all_path_dis = []
		title_change = config['title_change'] if 'title_change' in config else {"chrom": "#chr"}
		config_change = dict([(v.lower(), k.lower()) for k, v in title_change.iteritems()]) ##改变的位于顺序的第几个
		not_inter_function = ["utr-3","utr-5","no-change","promoter"]
		try:
			tgp_filter = float(config['autoInterp_freq_threshold']['TGP_AF'])
		except:
			tgp_filter = 0.01
		try:
			pvfd_filter = config['autoInterp_freq_threshold']['PVFD_AF']
		except:
			pvfd_filter = 0.01
		try:
			inter_function = config['inter_function']
		except:
			inter_function = ["cds-del", "cds-ins", "cds-indel", "cds-loss", "init-loss",
			                  "splice", "splice-3", "splice-5", "stop-gain", "stop-loss",
			                  "altstart", "frameshift", "knockout", "missense", "nonsense", "span"]
		try:
			filter_genes = config["filter_genes"]
		except:
			filter_genes = list()
		panel_genes = newborndb.genes['all']
		dis_genes = newborndb.genes['disease'] ##166个疾病基因 /ifs4/HST_5B/USERS/zhaocece/Production/newborn_datasets/disease_gene_v5
		drug_genes = newborndb.genes['drug']##药物基因drug_gene_v5
		enter_all = copy.copy(newborndb.entertainment)
		titles = dict()
		indb = dict()
		check_gene = defaultdict(set) #check_gene,基因与疾病组合键的包含的位点集合，那有没有可能是用药相关的位点都在里面
		for line in gzip.open(annofile):
			if line.startswith("##"):
				continue
			tmp_rows = line.strip().lower().split('\t')
			if not len(titles) and (title_change["chrom"] not in tmp_rows):
				continue
				##还没读取表头前的都舍弃
			elif title_change["chrom"] in tmp_rows: ##表头行
				for num in xrange(len(tmp_rows)):
					if tmp_rows[num] in config_change:
						titles[config_change[tmp_rows[num]]] = num
					else:
						titles[tmp_rows[num]] = num
			else:
				line = norm_xml_format(line.strip())##正常突变行，信息拆分记录
				rows = line.split('\t')
				try:
					gene_sym = rows[titles["genesymbol"]]
					if gene_sym == "NADKD1" : gene_sym = "NADK2"
					if gene_sym not in panel_genes and gene_sym != ".":
						continue
					if re.compile('h', re.I).match(rows[titles["zygosity"]]) and re.compile('X').search(rows[titles["chrom"]]) and gender_message(sample)=="M":
						rows[titles["zygosity"]] = "Hemi"
					gtype = rows[titles["zygosity"]]
					if re.compile('hom', re.I).match(gtype):
						gtype = "Hom"
					elif re.compile('het', re.I).match(gtype):
						gtype = "Het"
#					else:
#						continue

					chrom = 'chr' + re.sub('^chr', '', rows[titles["chrom"]])
					start = rows[titles["start"]]
					stop = rows[titles["stop"]]
					ref = rows[titles["ref"]]
					call = rows[titles["call"]]
					mut_function = rows[titles["functionname"]]
					rows[titles["flank"]] = newborndb.flank[chrom.lower()][stop]
					try:
						if "tgdaf" in titles:
							tgf = float(rows[titles["tgdaf"]])
						else:
							tgf = float(rows[titles["pvfdaf"]])
					except:
						tgf = 0.0
					try:
						if "pvfdaf" in titles:
							pvfd = float(rows[titles["pvfdaf"]])
						else:
							pvfd = tgf
					except:
						##点float不成功，则设置为0
						pvfd = 0.0
					hgmd = rows[titles['hgmdpred']].lower() if 'hgmdpred' in titles else "dm"
					try:
						trans = rows[titles["transcript"]]
					except:
						trans = "."
					try:
						chgvs = rows[titles["chgvs"]]
					except:
						chgvs = "."
					annotation_key = (chrom, start, stop, ref, call)
					enter_k1 = (chrom, start, stop, ref, call, "Ref")
					enter_k2 = (chrom, start, stop, ref, call, "Het")
					enter_k3 = (chrom, start, stop, ref, call, "Hom")
					enter_key = (chrom, start, stop, ref, call, gtype) ##gtype不就是Het或者Hom吗
					if  mut_function in not_inter_function and gene_sym not in drug_not_dis: continue  ##什么逻辑？在drug_and_dis里的不考虑的功能则会过滤掉，是的就是这个
					if chgvs and mut_function == 'intron' : ##intron保留10bp以内的

						chgvs_m = re.compile(r'c\.(\d+)[-\+](\d+)[A-Z]&gt;[A-Z]').search(chgvs)
						if chgvs_m:
							if int(chgvs_m.group(2)) > 10 and gene_sym not in drug_not_dis: continue

					if annotation_key in newborndb.mutdb:##如果突变是在已知数据库里的，则用已知数据库里的各种信息
						dbtrans = newborndb.mutdb[annotation_key][newborndb.database_title["nbtrans"]]
						gene_name = newborndb.mutdb[annotation_key][newborndb.database_title["#gene"]]
						disease = newborndb.mutdb[annotation_key][newborndb.database_title["disease"]].lower()
						if annotation_key not in indb or trans == dbtrans: ##indb是个临时汇总记录该样本在已有数据库的位点字典，如果没有记录过或者转录本是目的转录本，则记录,键是chr start end ref alt，感觉trans判定有点重复判定的意思，如果在数据库中并且也是这个转录本的不是已经记录过的
							indb[annotation_key] = rows
						if gtype == "Het":##为什么只有杂合才记录基因和疾病涉及的位点？
							tmp_key = (gene_name, disease)
							check_gene[tmp_key].add(annotation_key) ##check_gene基因与疾病组合键的包含的位点集合，这个变量的作用是啥？
					elif ((re.compile('dm', re.I).search(hgmd) or ((tgf <= tgp_filter) and (pvfd <= pvfd_filter)
					                                               and (mut_function in inter_function)))
					      and (gene_sym not in filter_genes) and (gene_sym in dis_genes)):
						##如果不在已知的数据库里，根据hgmd预测为有害或者频率很低，并且基因是疾病基因，不是需要过滤的基因，则往下走
						not_tmp = [gene_sym, trans, chgvs, chrom]
						for i, j in sorted(newborndb.database_title.iteritems(), key=lambda s: int(s[1]))[4:]:
							if i in titles:
								not_tmp.append(rows[titles[i]])
							else:
								not_tmp.append(".")
						annotations["not_in_db"].append(back_xml_format("\t".join(not_tmp)) + "\n") ##按照预先设置的格式存到不在数据库中的信息里

					if enter_key in newborndb.entertainment: #enter_key = (chrom, start, stop, ref, call, gtype)
						annotations["entertainment"].extend(enter_all[enter_key]) ### iwant zhidao 这个是用来干嘛的，春娜说是个体特征
						if enter_k1 in enter_all: del enter_all[enter_k1]
						if enter_k2 in enter_all: del enter_all[enter_k2]
						if enter_k3 in enter_all: del enter_all[enter_k3]
				except Exception, e:
					logger.info("%s is not well-formed! %s" % (sample, e))
					raise
		tmp_out = os.path.join(outdir, sample + "_mut_indb")
		annotations["in_db"], drug_usage = parsing_mut_indb(batch_num, sample, indb, check_gene, message['gender'], newborndb, haplotype,
		                                                    titles, tmp_out, all_indb_out)
		####here key important
		##batch_num 批次或者项目编号
		##sample 样品名称
		##indb是个临时汇总记录的字典，如果不在外部数据库或者转录本是目的转录本，则记录
		####check_gene基因与疾病组合键的包含的位点集合
		##message['gender'] 性别信息
		##newborndb是我也记不起来的数据库内容
		##haplotype 项目分析的result路径下有id.xls 以及每个样本的好像是基因型但是就两行
		##titles bedanno的表头，直接看起来就是，应该没有其他改动
		##tmp_out
		##all_indb_out在已有的数据库中的记录

		fout = open(os.path.join(outdir, sample + "_mut_Notindb"), 'w')##不在数据库中的位点输出到文件，但实际运行中最后又删除了，annotations["not_in_db"]是中间记录的变量，annotation记录在库不在库以及各种库的位点情况
		database_title = sorted(newborndb.database_title.keys(), key=lambda s: int(newborndb.database_title[s])) ###database_title过去突变位点的详细注释的表头的对应哪个列的记录
		contral_varsheet_dict = dict()
		for i, j in varsheet_title.iteritems():
#		for i, j in newborndb.database_title.iteritems():
			if i.lower() in varsheet_dict:
				tmp_weight = varsheet_dict[i.lower()]["width"] or None
				tmp_otherdict = {n:m for n, m in varsheet_dict[i.lower()].iteritems() if n != "width" and m != 0}
				if tmp_weight or len(tmp_otherdict):
					contral_varsheet_dict[int(j)] = (tmp_weight, tmp_otherdict)
		fout.write("Project" + '\t' + "Sample" + '\t' + back_xml_format("\t".join([newborndb.dbrawtitle[i] for i in database_title if newborndb.database_title[i] <= newborndb.database_title["phgvs"]])))
		fout.write('\t'+"ForPrimerDesign" + '\t' + back_xml_format("\t".join([newborndb.dbrawtitle[i] for i in database_title if newborndb.database_title[i] >newborndb.database_title["phgvs"]])) + '\n')
#		fout.write(back_xml_format("\t".join([newborndb.dbrawtitle[i] for i in database_title])) + '\n')
		for line in annotations["not_in_db"]:
			rows = line.strip().split('\t')
			gene_sym = rows[newborndb.database_title["#gene"]]
			trans = rows[newborndb.database_title["nbtrans"]]
			chgvs = rows[newborndb.database_title["nbchgvs"]]
			phgvs = rows[newborndb.database_title["phgvs"]]
			funcregion = rows[newborndb.database_title["funcregion"]]
			zygosity = rows[newborndb.database_title["zygosity"]]
			chrom=line[5].lower() if re.search('chr',line[5],re.I) else 'chr'+str(line[5]) #不知道为啥是line[5]
			flank = rows[newborndb.database_title["flank"]]
			forPrimerDesign = (gene_sym,trans,chgvs,phgvs,funcregion,zygosity,flank)
			fout.write(batch_num + '\t' + sample + '\t')
			all_Notindb_out.write(batch_num + '\t' + sample + '\t')
			for i in range(len(rows)):
				if(i == newborndb.database_title["phgvs"]):
					fout.write(back_xml_format(rows[i]) + '\t' + back_xml_format("; ".join(forPrimerDesign)) + '\t')
					all_Notindb_out.write(back_xml_format(rows[i]) + '\t' + back_xml_format("; ".join(forPrimerDesign)) + '\t')
				else:
					fout.write(back_xml_format(rows[i]) + '\t')
					all_Notindb_out.write(back_xml_format(rows[i]) + '\t')
			fout.write('\n')
			all_Notindb_out.write('\n')
#		fout.writelines(annotations["not_in_db"]) 增加sample的信息
		fout.close()

		last_reference = list()
		message_header = base.message.replace('$sample_id$',re.sub('-\d','',sample)).replace('$guardian_name$',message['guardian_name']).replace('$name$',message['name'])
		tmp_header = re.findall(r'>\$([^<>]+)\$<', message_header)
		for word in set(tmp_header):
			message_header = re.sub('\$'+word+'\$', message[word], message_header)
		rep.writelines(message_header)
		#rep.writelines(base.header)
		disease_definition = []
		genetic_suggestions = []
		dis_all = []
		references_disease_count = 0
		reference_dic = {}
		if len(annotations["in_db"]["pathgenic"])  or len(annotations["in_db"]["carrier"]) : ##weishenmene
			disease_count = set()
			dis_main_tmp = list()

			header = base.path_table_discripe.replace(r'1、', '')
			#header = header.replace(r'2、未检测出与遗传病相关的已知致病突变。检测结果提示未见有相关疾病患病风险和致病突变携带情况。',\
			#	'2、未检测出与遗传病相关的致病突变。检测结果提示未见有相关疾病患病风险和致病突变携带情况。')
			if len(annotations["in_db"]["pathgenic"]):
				dis_tmp = set()
				for rows in sorted(annotations["in_db"]["pathgenic"]):
					disease_ch, disease_en, gene_name, inheritance, chgvs, phgvs, functions, gtype, pathogenic, trans ,references = rows
					chgvs = get_alt(chgvs)
					if pathogenic == '已知致病':
						all_path_dis.append([disease_ch,'P','阳性'])
					else:
						all_path_dis.append([disease_ch,'LP','阳性'])
					try:
						last_table[partlast[disease_en.lower()][gene_name.lower()]] = True
					except:
						continue
					disease_count.add(disease_ch)
					if gtype == "纯合":
						gtype = "Hom"
					elif gtype == "杂合":
						gtype = "Het"

					tableline = ''.join([base.path_table_discr, "患病", base.disease_main_table_dis, newborndb.gene_diseasedb[gene_name][1], 
						base.disease_main_table_ini, newborndb.gene_diseasedb[gene_name][2], base.disease_main_table_gene, gene_name+'（'+trans+'）',
						base.disease_main_table_cds, chgvs,
						base.disease_main_table_phgvs, re.sub('p.0\?','.',phgvs),base.disease_main_table_function, functions,
						base.disease_main_table_gtype, gtype, base.disease_main_table_path,pathogenic, base.disease_main_table_end])
					dis_main_tmp.append(tableline)
					dis_all.append(newborndb.gene_diseasedb[gene_name][1])
					if newborndb.gene_diseasedb[gene_name][1] not in dis_all :
						dis_all.append(newborndb.gene_diseasedb[gene_name][1])

					temp_len = len(re.findall(r'\[\d+\]',references))
					if temp_len >= 1 :
						references_temp_list = re.split(r'\[\d+\]',references)
						del references_temp_list[0]
						print references_temp_list
						nums_temp = []
						for refer in references_temp_list:
							refer = re.sub('^\s+|\s+$','',refer)
							if refer in reference_dic:
								nums_temp.append(reference_dic[refer])
								continue
							elif len(refer) :
								references_disease_count = references_disease_count + 1 
								reference_dic[refer] = str(references_disease_count)
								nums_temp.append(reference_dic[refer])
								last_reference.append(refer)
								#references_pa_ca.append(refer)
						num = ','.join(nums_temp)
						if ('pathogenic',functions) in newborndb.definition_config:
							temp = newborndb.definition_config[('pathogenic',functions)][2].replace('$gene',gene_name).replace('$muntation$',chgvs).\
							replace('$disease$',disease_ch).replace('$type$',pathogenic).replace('$num$',num)
							temp = base.disease_main_definition_start+temp+base.disease_main_definition_end
							if not temp in disease_definition : 
								disease_definition.append(temp)
					else:
						if ('pathogenic',functions) in newborndb.definition_config:
							temp = newborndb.definition_config[('pathogenic',functions)][3].replace('$gene',gene_name).replace('$muntation$',chgvs).\
							replace('$disease$',disease_ch).replace('$type$',pathogenic)
							temp = base.disease_main_definition_start+temp+base.disease_main_definition_end
							if not temp in disease_definition : 
								disease_definition.append(temp)					
					dis_tmp.add(disease_ch)
				header = header.replace('$dis_of_path$','、'.join(dis_tmp))
			else:
				header = header.replace('$dis_of_path$患病风险和','')
			if len(annotations["in_db"]["carrier"]):
				risk_tmp = set()
				for rows in sorted(annotations["in_db"]["carrier"]):
					disease_ch, disease_en, gene_name, inheritance, chgvs, phgvs, functions, gtype, pathogenic, trans ,references = rows
					chgvs = get_alt(chgvs)
					disease_count.add(disease_ch)

					if gtype == "纯合":
						gtype = "Hom"
					elif gtype == "杂合":
						gtype = "Het"
					if pathogenic == '已知致病':
						all_path_dis.append([disease_ch,'P','携带'])
					else:
						all_path_dis.append([disease_ch,'LP','携带'])
					tableline = ''.join([base.path_table_discr, "携带者，对后代产生影响",base.disease_main_table_dis, newborndb.gene_diseasedb[gene_name][1], 
						base.disease_main_table_ini, newborndb.gene_diseasedb[gene_name][2], base.disease_main_table_gene, gene_name+'（'+trans+'）',
						base.disease_main_table_cds, chgvs,
						base.disease_main_table_phgvs, re.sub('p.0\?','.',phgvs),base.disease_main_table_function, functions,
						base.disease_main_table_gtype, gtype, base.disease_main_table_path,pathogenic, base.disease_main_table_end])
					dis_main_tmp.append(tableline)
					suggestion_xr = base.disease_main_definition_start+'携带者本人一般不会发病。但受检者'\
							'男性后代有1/2的风险为患者，女性后代有1/2的风险为携带者。建议受检者配偶婚孕前进行'\
							'相关疾病基因检测及遗传咨询。'+base.disease_main_definition_end
					suggestion_ar =base.disease_main_definition_start + '携带者本人一般不会发病。当携带者的未来的配偶也携带一样的突变时，他们的孩子将有25%的可'\
							'能患有该疾病，建议受检者配偶婚孕前进行相关疾病基因检测及遗传咨询。'+base.disease_main_definition_end
					
					if re.compile('AR').search(newborndb.gene_diseasedb[gene_name][2]) and suggestion_ar not in genetic_suggestions :
						genetic_suggestions.append(suggestion_ar)
					elif newborndb.gene_diseasedb[gene_name][2] == 'XR' and suggestion_xr not in genetic_suggestions :
						genetic_suggestions.append(suggestion_xr)

					temp_len = len(re.findall(r'\[\d+\]',references))
					if temp_len >= 1 :
						references_temp_list = re.split(r'\[\d+\]',references)
						del references_temp_list[0]
						print references_temp_list
						nums_temp = []
						for refer in references_temp_list:
							refer = re.sub('^\s+|\s+$','',refer)
							if refer in reference_dic:
								nums_temp.append(reference_dic[refer])
								continue
							elif len(refer) :
								references_disease_count = references_disease_count + 1 
								reference_dic[refer] = str(references_disease_count)
								nums_temp.append(reference_dic[refer])
								last_reference.append(refer)
								#references_pa_ca.append(refer)
						num = ','.join(nums_temp)
						if ('pathogenic',functions) in newborndb.definition_config:
							temp = newborndb.definition_config[('pathogenic',functions)][2].replace('$gene',gene_name).replace('$muntation$',chgvs).\
							replace('$disease$',disease_ch).replace('$type$',pathogenic).replace('$num$',num)
							if not temp in disease_definition :
								temp = base.disease_main_definition_start+temp+base.disease_main_definition_end 
								disease_definition.append(temp)
					else:
						if ('pathogenic',functions) in newborndb.definition_config:
							temp = newborndb.definition_config[('pathogenic',functions)][3].replace('$gene',gene_name).replace('$muntation$',chgvs).\
							replace('$disease$',disease_ch).replace('$type$',pathogenic)
							temp = base.disease_main_definition_start+temp+base.disease_main_definition_end
							if not temp in disease_definition : 
								disease_definition.append(temp)
					
					risk_tmp.add(disease_ch)
				header = header.replace('$dis_of_risk$','、'.join(risk_tmp))
			else:
				header = header.replace('和$dis_of_risk$的致病突变携带情况','')
			header = header.replace(r'$num_of_disease$',str(len(annotations["in_db"]["pathgenic"])+len(annotations["in_db"]["carrier"])))
			rep.writelines(header)
			rep.writelines(base.table_start)
			rep.writelines(base.path_table_title)
			rep.writelines("".join(dis_main_tmp))
			rep.writelines(base.table_end)
			
#			if len(annotations["in_db"]["carrier"]): rep.writelines(base.carrier)
		else:
			header = base.path_table_discripe.replace(r'1、共检测出$num_of_disease$个与遗传病相关的已知致病突变，检测结果提示'
			                                          r'有$dis_of_path$患病风险和$dis_of_risk$的致病突变携带情况。'
			                                          r'详见下表：2、','')
			rep.writelines(header)
			rep.writelines(base.add_blank)


		rep.writelines(base.disease_main_definition_title)
		rep.writelines(''.join(disease_definition))
		rep.writelines(base.disease_main_genetic_suggestions)
		
		if len(dis_all):
			dis = '、'.join(dis_all)
			genetic_suggestions_o = (base.disease_main_definition_start+'结果提示有%s患病风险，建议进行遗传咨询，'\
				'并进行临床相关检查。'+base.disease_main_definition_end)%dis
			rep.writelines(genetic_suggestions_o)
		
		if len(genetic_suggestions): rep.writelines(''.join(genetic_suggestions))
		rep.writelines(base.disease_main_background_title)


		for disease_ch, disease_en, mut_definiton, references in sorted(annotations["in_db"]["mut_definiton"]):
			if not re.compile(r'\[\d+\]').search(references):continue
			for termp_mut_refrence in re.split('\[[0-9]+\]', references):
					termp_mut_refrence = re.sub('^\s+|\s+$','',termp_mut_refrence)
					if termp_mut_refrence not in last_reference and len(termp_mut_refrence) > 3:
						last_reference.append(termp_mut_refrence)

		for disease_ch, disease_en, mut_message, references in sorted(annotations["in_db"]["disease_message"]):
			for termp_mut_refrence in re.split('\[[0-9]+\]', references):
				termp_mut_refrence = re.sub('^\s+|\s+$','',termp_mut_refrence)
				if termp_mut_refrence not in last_reference and len(termp_mut_refrence) > 3:
					last_reference.append(termp_mut_refrence)
			temp_mut_message = re.sub('\[[^A-Za-z\[\]]+\]','',mut_message).replace('$$',base.newline)
			rep.writelines(base.disease_main_background_start)
			rep.writelines(disease_ch + "（" + str(disease_en) + "）" )
			rep.writelines(base.disease_main_background_middle)
			rep.writelines(temp_mut_message)
			rep.writelines(base.disease_main_background_end)
		

		if len(annotations["in_db"]["drug"]): ##kaishikan
			drug = base.drug.replace(r'$num_of_drug$',str(len(drug_usage.keys()))).replace(r'1、', '')
			drug = drug.replace(r'2、在检测的药物基因组范围内，检测结果提示未见有相关药物剂量调整或换药情况。','')
			rep.writelines(drug)
			rep.writelines(base.table_start)
			rep.writelines(base.drug_table_title)
			for rows in sorted(annotations["in_db"]["drug"]):
				drug_ch, drug_en, drug_gene, haps, drug_description, drug_short_ch = rows
				try:
					last_table[partlast[drug_en.lower()][drug_gene.lower()]] = True
				except:
					continue
				if len(drug_usage[drug_en]) > 1: drug_short_ch = "药物剂量请参考国际标准"
				tableline = ''.join([base.drug_table_drug, drug_ch,base.drug_table_gene, drug_gene,base.drug_table_type,
				                     haps,base.drug_table_chtype, drug_description,base.drug_table_sug,
				                     drug_short_ch,base.drug_table_end])
				rep.writelines(tableline)
			rep.writelines(base.table_end)
		else:
			drug = base.drug.replace(r'1、检测结果提示有$num_of_drug$种药物需要调整剂量或换药，如果宝宝在以后的成长过程中，有需'
			                         r'要用到这$num_of_drug$种药物的情况出现，请您务必将检测结果提供给临床医生，以便为宝宝提供个体化用药。检测'
			                         r'结果详见下表：2、','')
			rep.writelines(drug)

		#total_path = copy.copy(annotations["in_db"]["pathgenic"])
		#total_path.update(annotations["in_db"]["carrier"])		
		rep.writelines(base.drug_main_definition_title)
		drug_definiton = list()
		for drug_ch, drug_en, drug_long_ch, references in sorted(annotations["in_db"]["drug_definiton"]):
			for termp_mut_refrence in re.split('\[[0-9]+\]', references):
				termp_mut_refrence = re.sub('^\s+|\s+$','',termp_mut_refrence)
				if termp_mut_refrence not in last_reference and len(termp_mut_refrence) > 3:
					last_reference.append(termp_mut_refrence)

				temp_drug_definiton = re.sub('\[[^A-Za-z\[\]]+\]','',drug_long_ch)
				if temp_drug_definiton not in drug_definiton:
					drug_definiton.append(temp_drug_definiton)
		for temp in drug_definiton:
			rep.writelines(base.drug_main_definition_start)
			temp = temp.replace('$$',base.newline)
			rep.writelines(temp)
			rep.writelines(base.drug_main_definition_end)
		rep.writelines(base.drug_main_background)
		for drug_ch, drug_en, drug_back_ch, references in sorted(annotations["in_db"]["drug_message"]):
			for termp_drug_refrence in re.split('\[[0-9]+\]', references):
				termp_drug_refrence = re.sub('^\s+|\s+$','',termp_drug_refrence)
				if termp_drug_refrence not in last_reference and len(termp_drug_refrence) > 3:
					last_reference.append(termp_drug_refrence)
			temp_drug_message = re.sub('\[[^A-Za-z\[\]]+\]','',drug_back_ch).replace('$$',base.newline)
			rep.writelines(base.drug_main_background_start)
			rep.writelines(drug_ch + "（" + str(drug_en) + "）" )
			rep.writelines(base.drug_main_background_middle)
			rep.writelines(temp_drug_message)
			rep.writelines(base.drug_main_background_end)


		result_table = base.result_table
		tmp_foot = re.findall(r'>\$([^<>]+)\$<', result_table)
		for word in set(tmp_foot):
			word_tmp = r'$'+word.lower()+ r'$'
			if word_tmp in last_table:
				word = '</w:rPr><w:t>' + r'$'+word+r'$'
				result_table = result_table.replace(word, r'<w:color w:val="FF0000"/></w:rPr><w:t>（＋）')
			else:
				word = r'$'+word+r'$'
				if word == r'$'+'result_seal_user_flag_23'+r'$' or word == r'$'+'result_seal_user_flag_42'+r'$' or word == r'$'+'result_seal_user_flag_12'+r'$' : continue
				result_table = result_table.replace(word, r'（－）')
		rep.writelines(result_table)
		#rep.writelines(base.reference_message)
		for num, ref in enumerate(last_reference, 1):
			rep.writelines(base.reference_start)
#			rep.writelines("["+str(num)+"] ")
			rep.writelines(ref)
			rep.writelines(base.reference_end)
		rep.writelines(base.science_read)
		rep.writelines(base.id_main)
		id_result = base.id_result
		id_results = re.findall(r'>\$([^<>]+)\$<', id_result)
		for word in set(id_results):
			word_tmp = r'$'+word.upper()+ r'$'
			if word in gene_id:
				id_result = id_result.replace(word_tmp,gene_id[word])
			else:
				id_result = id_result.replace(word_tmp,"-")
		rep.writelines(id_result)

		raw_entertainment_main = base.entertainment_main
		xls_s1  = open(os.path.join(outdir, str(sample) + '_Characteristic'), 'w')
		xls_s2  = open(os.path.join(outdir,str(sample) + '_Character_uncover'), 'w')
		xls1 = set()
		xls2 = set()
		uncover = set()
		entertainments_gtype = defaultdict(list)
		enter_summary = defaultdict(list)
		enter_fre = defaultdict(list)
		xls_s1.writelines("Num\tCharacteristic\tShort\tResult\tReference\n")
		xls_s2.writelines("Num\tCharacteristic\tShort\tResult\tReference\n")
		for chrom, start, end, ref, call, gtype in enter_all.keys():
			new_key = (chrom, start, end, ref, call, "Ref")
			if depth_from_bam:
				cover = depth_from_bam.count_coverage(str(chrom), int(start), int(end))
				tmpdepth = sum([int(k[-1]) for k in cover]) if len(cover) else 0
			else:
				tmpdepth = 4
			for value in enter_all[new_key]:
				if value not in annotations["entertainment"]:
					annotations["entertainment"].append(value)
					num, g_type, result, fre = value[0], value[-3], value[-2], value[-1]
					if not tmpdepth:
						uncover.add(num)
		if depth_from_bam:
			depth_from_bam.close()
		for value in sorted(annotations["entertainment"]):
			num, g_type, result, fre = value[0], value[-3], value[-2], value[-1]
			entertainments_gtype[num].append(g_type)
			enter_summary[num].append(float(result))
			enter_fre[num].append(float(fre))

		for raw_num in entertainments_gtype.keys():
			bins = round(float(sum(enter_summary[raw_num]))/len(enter_summary[raw_num]), 3)
			fre = max(0.001, reduce(lambda x,y:x*y, enter_fre[raw_num]))
			if bins > 1.0:
				from_description = newborndb.description(raw_num, 'x>1')
			elif bins == 1.0:
				from_description = newborndb.description(raw_num, 'x=1')
			else:
				from_description = newborndb.description(raw_num, 'x<1')
			from_description = norm_xml_format(from_description)
			number, factor, gender, characteristic, trait, short, results, reference = from_description.split('\t')
			results = results + r'，该基因型在中国人群中的频率约为' + float_to_per(fre)
			tmp = "\t".join([number, characteristic, short, results, reference]) + "\n"
			tmp = back_xml_format(tmp)
			xls1.add(tmp)
			if number in uncover:
				xls2.add(tmp)
			if gender.lower() == "m":
				gender = "男"
			elif gender.lower() == "f":
				gender = "女"
			for g_type_num in range(len(entertainments_gtype[raw_num])):
				gt_num = g_type_num + 1
				raw_entertainment_main = re.sub('\$'+raw_num+'_Gtype_'+str(gt_num)+'\$',
				                                str(entertainments_gtype[raw_num][g_type_num]), raw_entertainment_main)
				raw_entertainment_main = re.sub('\$'+raw_num+'_Fre_'+str(gt_num)+'\$',
				                                str(enter_fre[raw_num][g_type_num]), raw_entertainment_main)
			if (gender == 'U' or gender == message['gender']) and number not in uncover:
				raw_entertainment_main = re.sub('\$'+raw_num+'_Res_1\$', short, raw_entertainment_main)
				raw_entertainment_main = re.sub('\$'+raw_num+'_Res_2\$', results, raw_entertainment_main)
			elif number in uncover:
				raw_entertainment_main = re.sub('\$'+raw_num+'_Res_1\$', characteristic, raw_entertainment_main)
				raw_entertainment_main = re.sub('\$'+raw_num+'_Res_2\$', r'＊该表征未覆盖＊', raw_entertainment_main)
			else:
				raw_entertainment_main = re.sub('\$'+raw_num+'_Res_1\$', characteristic, raw_entertainment_main)
				ket_word = "＊该表征为%s性特征＊" % gender
				raw_entertainment_main = re.sub('\$'+raw_num+'_Res_2\$', ket_word, raw_entertainment_main)
		#rep.writelines(raw_entertainment_main)
		#note_message = base.note_message.replace(r'$datetime$', report_times)
		#rep.writelines(note_message)
		#rep.writelines(base.foot)
		raw_entertainment_main = raw_entertainment_main.replace(r'$datetime$', report_times)
#		raw_entertainment_main = raw_entertainment_main.replace(r'datetime', report_times)
		rep.writelines(raw_entertainment_main)
		rep.close()
		for tmp in sorted(xls1):
			xls_s1.writelines(tmp)
		for tmp in sorted(xls2):
			xls_s2.writelines(tmp)
		xls_s1.close()
		xls_s2.close()

		try:
			parsefile(os.path.join(outdir, sample, 'word/document.xml'))
		except Exception, e:
			logger.info("%s is NOT well-formed! %s" % (sample,e))
			raise

		text = "-".join('-'.join(i) for i in all_path_dis) if len(all_path_dis) else "阴性"
		outdocx = os.path.join(outdir, str(sample) + '.docx')
		compile_dir(os.path.join(outdir, sample), outdocx)
		report_date = time.strftime('%Y%m%d',time.localtime(time.time()))
		text = re.sub(r'/','or',text)
		report_name_check = re.sub('-\d','',str(sample)) + '_' + "DX0630" + '_' + "高通量" + '_' + str(report_date) + '-' + '安馨可新生儿基因检测报告(VIP)-' + str(message['gender']) + '-' + str(text) + '.docx'
		if len(report_name_check) >=256:
			os.rename(outdocx,os.path.join(outdir, re.sub('-\d','',str(sample)) + '_' + "DX0630" + '_' + "高通量" + '_' + str(report_date) + '-' + '安馨可新生儿基因检测报告(VIP)-' +
		                                str(message['gender']) + '-' + '.docx'))
		else:
			os.rename(outdocx, os.path.join(outdir, re.sub('-\d','',str(sample)) + '_' + "DX0630" + '_' + "高通量" + '_' + str(report_date) + '-' + '安馨可新生儿基因检测报告(VIP)-' +
		                                str(message['gender']) + '-' + str(text) + '.docx'))
		rep.close()
		test2xlsx(os.path.join(outdir, str(sample) + '_final_result'),
	          (os.path.join(outdir, str(sample) + '_mut_indb'), contral_varsheet_dict),
	          (os.path.join(outdir,str(sample) + '_mut_Notindb'), contral_varsheet_dict),
	          (os.path.join(outdir,str(sample) + '_Characteristic'), None),
	          (os.path.join(outdir,str(sample) + '_Character_uncover'), None))
		for files in attach:
			if os.path.isfile(files) and not os.path.exists(os.path.join(outdir, os.path.basename(files))):
				shutil.copy2(files, outdir)
		compile_dir(outdir, os.path.join(os.path.dirname(outdir),"Newborn_Test_report_of_" + str(sample) + "_VIP_" + report_type + ".zip"))
		logger.info("  |_ _ _ Finish writing reports of %s " % sample)
	except:
		rep.close()
		shutil.rmtree(os.path.join(outdir, sample), ignore_errors=True)
		raise


def main():
	usage = "Usage: %prog -i [sample.list] [options] "
	author = "Author: huangzhiwei@genomics.cn"
	description = __doc__.strip()
	parser = OptionParser(usage=usage, version="%prog 5.2", description=description, epilog=author)
	expect = OptionGroup(parser, 'Expected arguments', 'Caution: These parameters are necessary for %pro.')
	expect.add_option('-i', '--list', help='The Sample List (required) ', dest='sample_list') ##shurule
	parser.add_option_group(expect)
	optinal = OptionGroup(parser, 'Optional arguments', 'Caution: If you do not set these parameters in addition,'
	                                                    ' %pro will select the default value.')
	optinal.add_option('-o', '--outdir', help='The output dir [ default: `pwd` ] ', dest='outdir',
	                   default=os.getcwd() + '/')
	optinal.add_option('-p', '--hap', help='The haplotype result dir ', dest='haplotype', default=None)##shurule
	optinal.add_option('-b', '--bam', help='The bam file dir which contain all samples\'s final bam ,'
	                                       '[default:`[sample_list]/../alignment]/`', default=None, dest='bamdir')
	optinal.add_option('-m', '--message', help='The sample message FILE , format must set in config json file',
	                   default=None, dest='message_file')##病人就医简单资料
	optinal.add_option('-a', '--attach', help='The attachment FILE , if more than two files, join them with ":" ',
	                   dest='attachment', default=str())##粘贴的附件
	optinal.add_option('-c', '--config', help='The config json file, '
	                                          '[default: /ifs5/ST_TRANS_CARDIO/PUB/DATA/newborn_datasets/pro/config.json]',
	                   dest='base_config', default=os.path.dirname(os.path.realpath(__file__)) + '/config.json')
	optinal.add_option('-v', '--vcfanno', help='The dir which contain the samples\'s annotation results, '
	                                           '[default: [samplist]/../annotation/variation/]', dest='vcfanno',
	                   default=None)
	optinal.add_option('-d', '--database', help='The dir contain all databases, [default]', dest='databases',
	                   default=os.path.dirname(os.path.realpath(__file__)) + '/../')
	optinal.add_option('-t', '--type', help='The newborn report version (v5/v1), [ default: v5 ]',
	                   dest='Types', choices=['v5'], default='v5')
	optinal.add_option('--debug', help='debug model', dest='debug', default=False, action="store_true")
	parser.add_option_group(optinal)
	(options, args) = parser.parse_args()

	if not options.sample_list:
		parser.print_help()
		return 1
	sample_list = os.path.realpath(options.sample_list)
	if haps_check.haps_check(sample_list) is False :  ##看一下是不是样本之间是对应的，haps里的和samplelist，列数是否对，每个样本的31个基因里的特定组合是否存在，以及每个样本的xls文件是否存在，这个特定组合和用药是相关的吗，每个样本的xls是用来干嘛的
		exit(0)
	outdir = os.path.realpath(options.outdir)
	bamdir = os.path.realpath(options.bamdir) if options.bamdir else os.path.join(
		os.path.dirname(sample_list), '../alignment/')
	haplotype = options.haplotype if options.haplotype  else os.path.join(os.path.dirname(sample_list), 'result') ##单体型的路径，干嘛用啊是不是用药！！！！
	base_config = os.path.realpath(options.base_config)
	vcfanno = os.path.realpath(options.vcfanno) if options.vcfanno else os.path.join(os.path.dirname(sample_list),
	                                                                                 '../annotation/variation/')  ##bedanno样本注释数据路径
	message_file = options.message_file ##实际没写这个参数
	databases = os.path.realpath(options.databases) ##数据库必须有
	report_type = options.Types.strip() ##输出报告的版本
	attachment = set(options.attachment.split(":"))##添加的附件
	brochure = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'xml/个体特征检测报告之我的基因秘密.pdf')
	model = options.debug
	if not model:
		attachment.add(brochure) ##各种附件资料
	if brochure not in attachment: ##必须存在那个宣传pdf
		logger.info('Annex lost, it seems to be not normal, this explanation appears only as a warning!!')
		logger.info('I suggest that you should use \'-a %s \' option !!' % brochure)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	temp_directory_name = tempfile.mkdtemp()
	base_xml = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'xml', 'xml_base_vip_' + report_type + '.pyc') ##承载信息的格式
	demo = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'xml', 'VIP_' + report_type + '.docx') ##word格式模板

	logger.info("Loadding the report demo, config !!!")

	demo_code, config, base = load_demo_config(demo, base_config, base_xml, temp_directory_name)##word格式模板,各种表头以及疾病和基因的对应关系？，承载信息的格式xml
	config文件包含的表
	# "message_title_num": "Explain the each column's description set in sample_message_file, 0 base",
	# "title_change": "If the annotation files has modified the title which on the left , change the right side into the new",
	# "autoInterp_freq_threshold": "Frequency filtering threshold",
	# "inter_function": "Functions need to be rendered",
	# "filter_genes": "genes should be filter out ",
	# "default_fam_varSheet_format": "Final excel's table format",
	# "pathogenic_Definition": "If the database set a new definition (left), please select the relevant on the right",
	# "disease_drug_table": "Developer Mode"
	print base
	varsheet_dict = defaultdict(dict) ##存变异excel输出表头
	varsheet_title = defaultdict(dict)
	try:
		varSheet_formact = open(config["default_fam_varSheet_format"], "r")
		varSheet_lines = varSheet_formact.readlines()
		varSheet_title = varSheet_lines.pop(0).strip() if varSheet_lines[0].startswith("#title") else "#title\twidth\thidden\tlevel"
		vs_title = {j:i for i,j in enumerate(varSheet_title.split("\t"))}
		varSheet_num = 0
		for varline in varSheet_lines:
			vs_rows = varline.strip().split("\t")
			if len(vs_rows) != len(vs_title):
				continue
			k_w = re.sub("^\s+|\s+$", "", vs_rows[vs_title["#title"]].lower())
			varsheet_dict[k_w] = {i:int(vs_rows[j]) for i,j in vs_title.iteritems() if i != "#title"} ##根据default_fam_varSheet_format存储每个表头列宽等信息
			varsheet_title[k_w] = varSheet_num ##根据default_fam_varSheet_format存储每个表头在哪一列
			varSheet_num += 1
		varSheet_formact.close()
	except:
		pass
	logger.info("Parsing databases !")
	if 'pathogenic_Definition' in config:
		pathogenic_token = set(config['pathogenic_Definition'].keys())
	else:
		pathogenic_token = {"pathogenic", "dm", "likely pathogenic", "dp", "drug"}
	#newborndb = ParsingDatabase(databases, report_type, pathogenic_token)
	newborndb = ParsingDatabase(databases, report_type) ##解析各种数据库，包括基因或者突变的详细信息，这里面kanzhe 包含了药物数据库，这个important

	logger.info("Parsing sample message!")
	if 'message_title_num' in config:
		#message_title_num病人送检单上的各种信息
		#message_file没有输入

		sample_message = sample_messages(sample_list, message_file, config['message_title_num']) ##函数的作用是生成sample_message？
	else:
		sample_message = sample_messages(sample_list, message_file)

	logger.info("Loadding GeneID result !!!")
	project_name = os.path.basename(os.path.dirname(os.path.realpath(bamdir))) #newbornV5_sqx18SZ000043_BGISEQ-50018SZ0003477_20180607_22
	result_all = zipfile.ZipFile(os.path.join(outdir, str(project_name) + '_VIP_Report_' + report_type +'.zip'), 'w')##该批次的打包文件
	gene_id_file = os.path.join(haplotype, 'id') if os.path.isdir(haplotype) else None
	gene_id_dict = sample_gene_id(gene_id_file)##读取id.xls的内容，记录成每个样本在每个基因（或者基因组合）的基因型？比如g-T-A g-C-G，并返回
	logger.info("Running, bugs please mail to huangyushan@genomics.cn")
	result_list = list()
	if "_" in project_name:
		project = project_name.strip().split("_")
		batch_num = project[1] if project[1] else "NewBorn" + report_type
	else:
		batch_num = "NewBorn" + report_type
	all_mut_indb = os.path.join(outdir, "samples_mut_indb")  ##为啥没这俩文件呢，没搜到
	all_mut_Notindb = os.path.join(outdir, "samples_mut_Notindb")
	if os.path.isfile(all_mut_indb):
		os.remove(all_mut_indb)
		all_indb_out = open(all_mut_indb, 'a') 
	else:
		all_indb_out = open(all_mut_indb, 'a')
	if os.path.isfile(all_mut_Notindb):
		os.remove(all_mut_Notindb)
		all_Notindb_out = open(all_mut_Notindb, 'a')
	else:
		all_Notindb_out = open(all_mut_Notindb, 'a')

	back_xml_format = NormXmlFormat(False)##是否标准化格式
	database_title = sorted(newborndb.database_title.keys(), key=lambda s: int(newborndb.database_title[s]))###过去的记录，已经解读过的位点的详细信息的表头，键是表头，值是表头所在的列，目的是将现在新样本的点以过去数据库的格式存储起来吧，所以才这么绕，在不在已有数据库里的存储格式都是一样的，在这是把表头先写好
	all_indb_out.write("Project" + '\t' "Sample" + '\t' + back_xml_format("\t".join([newborndb.dbrawtitle[i] for i in database_title if newborndb.database_title[i] <= newborndb.database_title["phgvs"]]))) ##phgvs
	all_indb_out.write('\t'+"ForPrimerDesign" + '\t' + back_xml_format("\t".join([newborndb.dbrawtitle[i] for i in database_title if newborndb.database_title[i] >newborndb.database_title["phgvs"]])) + '\n')
	all_Notindb_out.write("Project" + '\t' "Sample" + '\t' + back_xml_format("\t".join([newborndb.dbrawtitle[i] for i in database_title if newborndb.database_title[i] <= newborndb.database_title["phgvs"]])))
	all_Notindb_out.write('\t'+"ForPrimerDesign" + '\t' + back_xml_format("\t".join([newborndb.dbrawtitle[i] for i in database_title if newborndb.database_title[i] >newborndb.database_title["phgvs"]])) + '\n')

#	pool = multiprocessing.Pool(processes=min(cpu_count(), len(sample_message.keys())))
	for sample in sorted(sample_message.keys()):
		sample_outdir = os.path.join(outdir, sample)
		if not os.path.exists(sample_outdir):
			os.makedirs(sample_outdir)
		sample_hla = os.path.join(haplotype, sample + '_type') if os.path.isdir(haplotype) else None
#		pool.apply_async(auto_report, (sample, vcfanno, bamdir, temp_directory_name, sample_outdir, newborndb,
#		                               sample_hla, sample_message[sample], gene_id_dict[sample], config, base, attachment))
		auto_report(batch_num, sample, vcfanno, bamdir, temp_directory_name, sample_outdir, newborndb, sample_hla,
		                                      sample_message[sample], gene_id_dict[sample], config, base, attachment, 
											  report_type, varsheet_dict, varsheet_title, all_indb_out, all_Notindb_out)
		##auto_report用来生成每个样本的各种报告，各种所需的文件，并打包，生成格式化的vcf注释文件
		##base是xml的承载信息变量
		file_out = os.path.join(outdir,"Newborn_Test_report_of_" + str(sample) + "_VIP_" + report_type + ".zip")
		result_list.append(file_out)
	all_indb_out.close()
	all_Notindb_out.close()
	contral_varsheet_dict = dict()
	for i, j in varsheet_title.iteritems():##根据default_fam_varSheet_format存储每个表头在哪一列
		if i.lower() in varsheet_dict:##根据default_fam_varSheet_format存储每个表头列宽等信息
			tmp_weight = varsheet_dict[i.lower()]["width"] or None
			tmp_otherdict = {n:m for n, m in varsheet_dict[i.lower()].iteritems() if n != "width" and m != 0}
			if tmp_weight or len(tmp_otherdict):
				contral_varsheet_dict[int(j)] = (tmp_weight, tmp_otherdict)
	test2xlsx(os.path.join(outdir, 'all_samples_final_result'),
			(os.path.join(outdir, 'samples_mut_indb'), contral_varsheet_dict),
			(os.path.join(outdir, 'samples_mut_Notindb'), contral_varsheet_dict))
	result_list.append(os.path.join(outdir, 'all_samples_final_result.xlsx'))
	logger.info("Clean temp files !")
#	pool.close()
#	pool.join()
	for files in result_list:
		abs_path = os.path.realpath(files)
		rel_path = os.path.relpath(abs_path,os.path.dirname(outdir + '/'))
		result_all.write(abs_path, rel_path, zipfile.ZIP_STORED)
		os.remove(abs_path)
	result_all.close()

	logger.info("Jobs Finished, results stored in %s !" % os.path.join(outdir, str(project_name) + '_VIP_Report_' + report_type + '.zip'))
	return 0
if __name__ == '__main__':
	sys.exit(logger.info("Process finished with exit code : %d" % main()))
