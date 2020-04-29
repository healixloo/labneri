# ------------------------------------------------------------
# Author: Natalie Romanov
# ------------------------------------------------------------

class step1_string_preparation:

	@staticmethod
	def execute(**kwargs):
		folder = kwargs.get('folder','PATH')		
		species = kwargs.get('species','human')

		print('STEP1: load_string_data')
		step1_string_preparation.load_string_data(species, folder)
		print('STEP1: map_string_data_to_symbols')
		step1_string_preparation.map_string_data_to_symbols(stringData, folder, species)
		print('STEP1: prepare_string_json_dictionaries')
		step1_string_preparation.prepare_string_json_dictionaries(stringData, folder)

	@staticmethod
	def load_string_data(species, folder):
		#downloaded from STRING database
		#---------------------------------------------------------------------
		if species == 'human':
			filename = 'stringData_9606.protein.links.detailed.v10.txt.gz'
		elif species == 'mouse':
			filename = 'stringData_10090.protein.links.detailed.v10.txt.gz'

		stringData = DataFrameAnalyzer.open_in_chunks(folder, filename, delim = ' ', header = None)
		stringData.columns = ['protein2','neighborhood','fusion','cooccurence',
							  'coexpression','experimental','database','textmining',
							  'combined_score']
		return stringData

	@staticmethod
	def map_string_data_to_symbols(stringData, folder, species):
		proteinList1  =  [str(item).split('.')[1] for item in filter(lambda a:str(a)! = 'nan',list(set(stringData.index)))]
		proteinList2  =  [str(item).split('.')[1] for item in filter(lambda a:str(a)! = 'nan',list(set(stringData.protein2)))]
		protein_list  =  list(set(proteinList1 + proteinList2))

		m  =  Mapper(protein_list,
					 input = 'string,ensembl.protein',
					 output = 'symbol',
					 species = species)

		hitDict = m.hitDict['symbol']

		pro_list1 = list(stringData.index)
		pro_list2 = list(stringData.protein2)
		sym_list1 = list()
		sym_list2 = list()
		count  =  0
		for pro1,pro2 in zip(pro_list1, pro_list2):
			if str(pro1)! = 'nan' and str(pro2)! = 'nan':
				pro1 = str(pro1).split('.')[1]
				pro2 = str(pro2).split('.')[1]
				get_protein1 = hitDict.get(pro1,0)
				get_protein2 = hitDict.get(pro2,0)
				if get_protein1!=0:
					sym_list1.append(hitDict[pro1])
				else:
					sym_list1.append('')
				if get_protein2!=0:
					sym_list2.append(hitDict[pro2])
				else:
					sym_list2.append('')
			else:
				sym_list1.append('')
				sym_list2.append('')
			count+= 1
		stringData['symbol1'] = pd.Series(sym_list1, index = stringData.index)
		stringData['symbol2'] = pd.Series(sym_list2, index = stringData.index)
		stringData.to_csv(folder + species.upper() + '_protein.links.detailed.v10.5_withSymbol.tsv.gz',
						  compression = 'gzip', sep = '\t')
		return stringData

	@staticmethod
	def prepare_string_json_dictionaries(stringData, folder):
		string_data = stringData.copy()
		string_data['protein1'] = string_data.index
		string_data['ensp_interaction'] = string_data['protein1'] + ':' + string_data['protein2']
		string_data['interaction'] = string_data['symbol1'] + ':' + string_data['symbol2']

		string_data1 = string_data.copy()
		string_data1.index = string_data.symbol1
		string_data1_500 = string_data1[string_data1['combined_score']> = 500]
		string_data1_700 = string_data1[string_data1['combined_score']> = 700]
		string_data1_exp500 = string_data1[string_data1['experimental']> = 500]
		string_data1_exp700 = string_data1[string_data1['experimental']> = 700]

		string_dict1_all = {k: list(v) for k,v in string_data1.groupby("symbol1")["symbol2"]}
		string_dict1_500 = {k: list(v) for k,v in string_data1_500.groupby("symbol1")["symbol2"]}
		string_dict1_700 = {k: list(v) for k,v in string_data1_700.groupby("symbol1")["symbol2"]}
		string_dict1_exp500 = {k: list(v) for k,v in string_data1_exp500.groupby("symbol1")["symbol2"]}
		string_dict1_exp700 = {k: list(v) for k,v in string_data1_exp700.groupby("symbol1")["symbol2"]}

		string_data2 = string_data.copy()
		string_data2.index = string_data.symbol2
		string_data2_500 = string_data2[string_data2['combined_score']> = 500]
		string_data2_700 = string_data2[string_data2['combined_score']> = 700]
		string_data2_exp500 = string_data2[string_data2['experimental']> = 500]
		string_data2_exp700 = string_data2[string_data2['experimental']> = 700]

		string_dict2_all = {k: list(v) for k,v in string_data2.groupby("symbol2")["symbol1"]}
		string_dict2_500 = {k: list(v) for k,v in string_data2_500.groupby("symbol2")["symbol1"]}
		string_dict2_700 = {k: list(v) for k,v in string_data2_700.groupby("symbol2")["symbol1"]}
		string_dict2_exp500 = {k: list(v) for k,v in string_data2_exp500.groupby("symbol2")["symbol1"]}
		string_dict2_exp700 = {k: list(v) for k,v in string_data2_exp700.groupby("symbol2")["symbol1"]}

		symbolList1 = list(string_data['symbol1'])
		symbolList2 = list(string_data['symbol2'])
		symbol_list = list(set(symbolList1 + symbolList2))

		string_dict_all = dict()
		string_dict_500 = dict()
		string_dict_700 = dict()
		string_dict_exp500 = dict()
		string_dict_exp700 = dict()
		for scount,sym in enumerate(symbol_list):
			get_sym1 = string_dict1_all.get(sym, 0)
			get_sym2 = string_dict2_all.get(sym, 0)
			get_sym1_500 = string_dict1_500.get(sym, 0)
			get_sym2_500 = string_dict2_500.get(sym, 0)
			get_sym1_700 = string_dict1_700.get(sym, 0)
			get_sym2_700 = string_dict2_700.get(sym, 0)
			get_sym1_exp500 = string_dict1_exp500.get(sym, 0)
			get_sym2_exp500 = string_dict2_exp500.get(sym, 0)
			get_sym1_exp700 = string_dict1_exp700.get(sym, 0)
			get_sym2_exp700 = string_dict2_exp700.get(sym, 0)

			protein_list = list()
			protein_list_500 = list()
			protein_list_700 = list()
			protein_list_exp500 = list()
			protein_list_exp700 = list()
			if get_sym1!= 0:
				protein_list.append(string_dict1_all[sym])
			if get_sym2!= 0:
				protein_list.append(string_dict2_all[sym])
			if get_sym1_500!= 0:
				protein_list_500.append(string_dict1_500[sym])
			if get_sym2_500!= 0:
				protein_list_500.append(string_dict2_500[sym])
			if get_sym1_700!= 0:
				protein_list_700.append(string_dict1_700[sym])
			if get_sym2_700!= 0:
				protein_list_700.append(string_dict2_700[sym])
			if get_sym1_exp500!= 0:
				protein_list_exp500.append(string_dict1_exp500[sym])
			if get_sym2_exp500!= 0:
				protein_list_exp500.append(string_dict2_exp500[sym])
			if get_sym1_exp700!= 0:
				protein_list_exp700.append(string_dict1_exp700[sym])
			if get_sym2_exp700!= 0:
				protein_list_exp700.append(string_dict2_exp700[sym])
			string_dict_all[sym] = list(set(utilsFacade.flatten(protein_list)))
			string_dict_500[sym] = list(set(utilsFacade.flatten(protein_list_500)))
			string_dict_700[sym] = list(set(utilsFacade.flatten(protein_list_700)))
			string_dict_exp500[sym] = list(set(utilsFacade.flatten(protein_list_exp500)))
			string_dict_exp700[sym] = list(set(utilsFacade.flatten(protein_list_exp700)))

		json.dump(string_dict_all, open(folder + 'string_interactors_geneName_ALL.json','wb'))
		json.dump(string_dict_500, open(folder + 'string_interactors_geneName_500.json','wb'))
		json.dump(string_dict_700, open(folder + 'string_interactors_geneName_700.json','wb'))
		json.dump(string_dict_exp500, open(folder + 'string_interactors_geneName_exp500.json','wb'))
		json.dump(string_dict_exp700, open(folder + 'string_interactors_geneName_exp700.json','wb'))

class tcga_ovarian_cancer_formatting_prep(object):
	def __init__(self, **kwargs):
		folder = kwargs.get('folder','PATH')

		print("load_data")
		self.data1,self.data2 = self.load_data1()

		print("merge_data")
		self.df = self.merge_data()
		print("load_metadata")
		self.metadata = self.load_metadata()

		print("load_ages_samples")
		self.ages,self.samples = self.load_ages_samples()
		
		print("load_new_data")
		self.df = self.load_data2()

		print("export_normalized_dataset")
		self.export_normalized_dataset()

		print("translate_data_ovarian")
		self.trans_dict = self.translate_data_ovarian()

		print("map_information_ovarian")
		self.data_ovarian = self.map_information_ovarian()

	def get_columns(self,data):
		cols = list()
		quant_cols = list()
		for col in list(data.columns):
			if col.find("Unshared") == -1:
				if col.find("Ratio")!= -1 and col.find("CONTROL") == -1:
					quant_cols.append(col)
				cols.append(col)
		return cols,quant_cols

	def load_data1(self):
		folder = self.folder
		fileName = "TCGA_Ovarian_PNNL_Proteome.itraq.tsv"
		data1 = DataFrameAnalyzer.getFile(folder,fileName)
		fileName = "TCGA_Ovarian_JHU_Proteome.itraq.tsv"
		data2 = DataFrameAnalyzer.getFile(folder,fileName)
		return data1,data2

	def load_data2(self):
		folder = self.folder
		ages,samples = self.ages,self.samples

		df = DataFrameAnalyzer.getFile(folder,"TCGA_Ovarian_ALL_Proteome.itraq.tsv")
		df = df.replace(np.nan,0.0)
		cols,qcols = self.get_columns(df)
		sample_cols = list()
		for sam in samples:
			sample = sam.split("TCGA-")[1]
			for col in qcols:
				if col.find(sample)!= -1:
					sample_cols.append(col)
		df = df.drop(["Mean","Median","StdDev"],0)
		quant_df = self.quantileNormalize(df)
		quant_df = quant_df.replace(0.0,np.nan)
		df = quant_df
		df = self.median_normalization(df)
		df = df[sample_cols]
		return df

	def merge_data(self):
		folder = self.folder
		data1,data2 = self.data1,self.data2

		cols1,qcols1 = self.get_columns(data1)
		cols2,qcols2 = self.get_columns(data2)
		qcols = list(set(qcols1).union(set(qcols2)))
		quant_data1 = data1[qcols1]
		quant_data1 = quant_data1.drop(["Mean","Median","StdDev"],0)
		qdict1 = quant_data1.to_dict()
		quant_data2 = data2[qcols2]
		quant_data2 = quant_data2.drop(["Mean","Median","StdDev"],0)
		qdict2 = quant_data2.to_dict()

		proteinList = list(set(data1.index).union(set(data2.index)))
		dfList = list()
		for protein in proteinList:
			tempList = list()
			for col in qcols:
				get_key1 = qdict1.get(col,"bla")
				get_key2 = qdict2.get(col,"bla")
				lst = list()
				if get_key1! = "bla" and get_key2! = "bla":
					get_protein1 = qdict1[col].get(protein,"bla")
					get_protein2 = qdict2[col].get(protein,"bla")
					if get_protein1! = "bla" and get_protein2! = "bla":
						tempList.append(np.mean([qdict1[col][protein],qdict2[col][protein]]))
					elif get_protein1! = "bla" and get_protein2 =  = "bla":
						tempList.append(qdict1[col][protein])
					elif get_protein1 =  = "bla" and get_protein2! = "bla":
						tempList.append(qdict2[col][protein])
					else:
						tempList.append("")
				elif get_key1! = "bla" and get_key2 =  = "bla":
					get_protein1 = qdict1[col].get(protein,"bla")
					if get_protein1! = "bla":
						tempList.append(qdict1[col][protein])
					else:
						tempList.append("")
				elif get_key1 =  = "bla" and get_key2! = "bla":
					get_protein2 = qdict2[col].get(protein,"bla")
					if get_protein2! = "bla":
						tempList.append(qdict2[col][protein])
					else:
						tempList.append("")
			dfList.append(tempList)
		df = pd.DataFrame(dfList)
		df.index = proteinList
		df.columns = qcols
		df.to_csv(folder + "TCGA_Ovarian_ALL_Proteome.itraq.tsv",sep = "\t")
		return df

	def quantileNormalize(self,df_input):
	    df = df_input.copy()
	    #compute rank
	    dic = {}
	    for col in df:
	        dic.update({col : sorted(df[col])})
	    sorted_df = pd.DataFrame(dic)
	    rank = sorted_df.mean(axis = 1).tolist()
	    #sort
	    for col in df:
	        t = np.searchsorted(np.sort(df[col]), df[col])
	        df[col] = [rank[i] for i in t]
	    return df

	def median_normalization(self,df_input):
		df = df_input.copy()
		median_list = list()
		for col in df:
			median_list.append(np.median(finite(df[col])))
		dic = dict()
		for c,col in enumerate(df.columns):
			dic.setdefault(col,[])
			dic[col] = list(np.array(df[col]) - median_list[c])
		med_df = pd.DataFrame(dic)
		med_df.index = df.index
		return med_df

	def load_metadata(self):
		folder = self.folder
		metadata = DataFrameAnalyzer.getFile(folder,"tcga_ovarian_metadata_patients.txt")
		ageList = list(metadata["days_to_birth"])
		new_ageList = list()
		new_diagnosisList = list()
		for days in ageList:
			try:
				age = float(abs(int(days)))/float(365.0)
				new_ageList.append(age)
			except:
				new_ageList.append(np.nan)
		metadata["age"] = pd.Series(new_ageList,index = metadata.index)
		metadata.to_csv(folder + "tcga_ovarian_metadata_patients.txt", sep = "\t")
		return metadata

	def load_ages_samples(self):
		metadata = self.metadata
		metadata.index = metadata["bcr_patient_barcode"]
		meta_age_dict = metadata["age"].to_dict()
		samples = list()
		ages = list()
		for k in meta_age_dict.keys():
			age = meta_age_dict[k]
			if str(age)! = "nan":
				ages.append(age)
			else:
				ages.append(0)
			samples.append(k)
		ages,samples = zip(*sorted(zip(ages,samples)))
		return ages,samples

	def export_normalized_dataset(self):
		df = self.df
		folder = self.folder
		df.to_csv(folder + "TCGA_Ovarian_ALL_Proteome.itraq_NORM.tsv",sep = "\t")

	def translate_data_ovarian(self):
		data_ovarian = self.df
		protein_list = list(data_ovarian.index)
		trans_df = mg.querymany(protein_list,
								scopes = 'entrezgene,ensembl.gene,symbol,reporter,accession,uniprot,name,alias,summary',
								fields = 'symbol,entrezgene,uniprot,summary,alias,go,pathway,name,homologene,homologene.id,ensembl.gene,ensembl.protein',
								species = "human",as_dataframe = True,returnall = True)
		trans_dict = trans_df["out"].to_dict()
		return trans_dict

	def map_information_ovarian(self):
		folder = self.folder
		data_ovarian = self.df
		trans_dict = self.trans_dict
		
		protein_list = list(data_ovarian.index)
		entrezList = list()
		ensembl_gene_list = list()
		ensembl_protein_list = list()
		nameList = list()
		aliasList = list()
		pathway_idList = list()
		pathway_nameList = list()
		uniprotList = list()
		for protein in protein_list:
			entrezList.append(trans_dict["entrezgene"][protein])
			nameList.append(trans_dict["name"][protein])
			if type(trans_dict["alias"][protein]) == list:
				aliasList.append(",".join(trans_dict["alias"][protein]))
			else:
				aliasList.append(trans_dict["alias"][protein])
			if str(trans_dict["uniprot"][protein])!= "nan":
				get_swissprot = trans_dict["uniprot"][protein].get("Swiss-Prot",0)
				get_trembl = trans_dict["uniprot"][protein].get("TrEMBL",0)
				uniList = list()
				if get_swissprot!= 0:
					if type(trans_dict["uniprot"][protein]["Swiss-Prot"]) == list:
						for item in trans_dict["uniprot"][protein]["Swiss-Prot"]:
							uniList.append(item)
					else:
						uniList.append(trans_dict["uniprot"][protein]["Swiss-Prot"])
				if get_trembl!= 0:
					if type(trans_dict["uniprot"][protein]["TrEMBL"]) == list:
						for item in trans_dict["uniprot"][protein]["TrEMBL"]:
							uniList.append(item)
					else:
						uniList.append(trans_dict["uniprot"][protein]["TrEMBL"])
				uniprotList.append(",".join(uniList))
			else:
				uniprotList.append("")
			if str(trans_dict["pathway"][protein])!= "nan":
				get_pathway = trans_dict["pathway"][protein].get("reactome",0)
				id_patList = list()
				name_patList = list()
				if get_pathway!= 0:
					if type(trans_dict["pathway"][protein]["reactome"]) == list:
						for item in trans_dict["pathway"][protein]["reactome"]:
							id_patList.append(item["id"])
							name_patList.append(item["name"])
					else:
						id_patList.append(trans_dict["pathway"][protein]["reactome"]["id"])
						name_patList.append(trans_dict["pathway"][protein]["reactome"]["name"])
				pathway_idList.append("//".join(id_patList))
				pathway_nameList.append("//".join(name_patList))
			else:
				pathway_idList.append("")
				pathway_nameList.append("")
			if str(trans_dict["ensembl"][protein])!= "nan":
				e_genes = list()
				e_proteins = list()
				try:
					get_gene = trans_dict["ensembl"][protein].get("gene",0)
					get_protein = trans_dict["ensembl"][protein].get("protein",0)
					if get_gene!= 0:
						if type(get_gene) == list:
							for item in get_gene:
								e_genes.append(item)
						else:
							e_genes.append(get_gene)
					if get_protein!= 0:
						if type(get_protein) == list:
							for item in get_protein:
								e_proteins.append(item)
						else:
							e_proteins.append(get_protein)
				except:
					get_gene = trans_dict["ensembl"][protein][0].get("gene",0)
					get_protein = trans_dict["ensembl"][protein][0].get("protein",0)
					if get_gene!= 0:
						if type(get_gene) == list:
							for item in get_gene:
								e_genes.append(item)
						else:
							e_genes.append(get_gene)
					if get_protein!= 0:
						if type(get_protein) == list:
							for item in get_protein:
								e_proteins.append(item)
						else:
							e_proteins.append(get_protein)
				ensembl_gene_list.append(",".join(e_genes))
				ensembl_protein_list.append(",".join(e_proteins))	
			else:
				ensembl_gene_list.append("")		
				ensembl_protein_list.append("")

		data_ovarian["entrez"] = pd.Series(entrezList,index = data_ovarian.index)
		data_ovarian["name"] = pd.Series(nameList,index = data_ovarian.index)
		data_ovarian["alias"] = pd.Series(aliasList,index = data_ovarian.index)
		data_ovarian["uniprot"] = pd.Series(uniprotList,index = data_ovarian.index)
		data_ovarian["pathway_id"] = pd.Series(pathway_idList,index = data_ovarian.index)
		data_ovarian["pathway_name"] = pd.Series(pathway_nameList,index = data_ovarian.index)
		data_ovarian["ensembl.gene"] = pd.Series(ensembl_gene_list,index = data_ovarian.index)
		data_ovarian["ensembl.protein"] = pd.Series(ensembl_protein_list,index = data_ovarian.index)

		proteinList = list(data_ovarian.index)
		complex_idList = list()
		complex_nameList = list()
		for protein in proteinList:
			complex_ids = list()
			complex_names = list()
			for complexID in complexDict:
				humanSubunits = complexDict[complexID]["humanGeneNames"]
				if protein in humanSubunits:
					complex_ids.append(complexID)
					complex_names.append(complexDict[complexID]["altName"][0])
			complex_idList.append(",".join(complex_ids))
			complex_nameList.append(",".join(complex_names))
		data_ovarian["complex_formal_ID"] = pd.Series(complex_idList,index = data_ovarian.index)
		data_ovarian["complexName"] = pd.Series(complex_nameList,index = data_ovarian.index)
		fileName = "TCGA_Ovarian_ALL_Proteome.itraq_NORM_INFO.tsv"
		data_ovarian.to_csv(folder + fileName, sep = "\t")
		return data_ovarian

class tcga_breast_cancer_formatting_prep(object):
	def __init__(self, **kwargs):
		folder = kwargs.get('folder','PATH')

		print("load_data")
		self.data = self.load_data1()

		print("merge_data")
		self.df = self.merge_data()

		print("load_metadata")
		self.metadata = self.load_metadata()

		print("load_ages_samples")
		self.ages,self.samples = self.load_ages_samples()
		
		print("load_new_data")
		self.df = self.load_data2()

		print("export_normalized_dataset")
		self.export_normalized_dataset()

		print("translate_data_breast")
		self.trans_dict = self.translate_data_breast()

		print("map_information_breast")
		self.data_ovarian = self.map_information_breast()

	def get_columns(self,data):
		cols = list()
		quant_cols = list()
		for col in list(data.columns):
			if col.find("Unshared") == -1:
				if col.find("Ratio")!= -1 and col.find("CONTROL") == -1:
					quant_cols.append(col)
				cols.append(col)
		return cols,quant_cols

	def load_data1(self):
		folder = self.folder
		fileName = "TCGA_Breast_BI_Proteome.itraq.tsv"
		data = DataFrameAnalyzer.getFile(folder,fileName)
		return data

	def merge_data(self):
		folder = self.folder
		data = self.data

		cols,qcols = self.get_columns(data)
		quant_data = data[qcols]
		quant_data = quant_data.drop(["Mean","Median","StdDev"],0)
		df = quant_data
		df.to_csv(folder + "TCGA_Breast_ALL_Proteome.itraq.tsv",sep = "\t")
		return df

	def quantileNormalize(self,df_input):
	    df = df_input.copy()
	    #compute rank
	    dic = {}
	    for col in df:
	        dic.update({col : sorted(df[col])})
	    sorted_df = pd.DataFrame(dic)
	    rank = sorted_df.mean(axis  =  1).tolist()
	    #sort
	    for col in df:
	        t = np.searchsorted(np.sort(df[col]), df[col])
	        df[col] = [rank[i] for i in t]
	    return df

	def median_normalization(self,df_input):
		df = df_input.copy()
		median_list = list()
		for col in df:
			median_list.append(np.median(finite(df[col])))
		dic = dict()
		for c,col in enumerate(df.columns):
			dic.setdefault(col,[])
			dic[col] = list(np.array(df[col]) - median_list[c])
		med_df = pd.DataFrame(dic)
		med_df.index = df.index
		return med_df

	def load_ages_samples(self):
		metadata = self.metadata
		meta_age_dict = metadata["age"].to_dict()
		samples = list()
		ages = list()
		for k in meta_age_dict.keys():
			age = meta_age_dict[k]
			if str(age)! = "nan":
				ages.append(age)
			else:
				ages.append(0)
			samples.append(k)
		ages,samples = zip(*sorted(zip(ages,samples)))
		return ages,samples

	def load_metadata(self):
		folder = self.folder
		metadata = DataFrameAnalyzer.getFile(folder, "tcga_brca_metadata_patients.txt")
		return metadata

	def load_data2(self):
		folder = self.folder
		ages,samples = self.ages,self.samples

		df = DataFrameAnalyzer.getFile(folder,"TCGA_Breast_ALL_Proteome.itraq.tsv")
		df = df.replace(np.nan,0.0)
		cols,qcols = self.get_columns(df)
		sample_cols = list()
		for sam in samples:
			sample = sam.split("TCGA-")[1]
			for col in qcols:
				if col.find(sample)!= -1:
					sample_cols.append(col)
		#df = df.drop(["Mean","Median","StdDev"],0)
		quant_df = self.quantileNormalize(df)
		quant_df = quant_df.replace(0.0,np.nan)
		df = quant_df
		df = self.median_normalization(df)
		df = df[sample_cols]
		return df
	
	def export_normalized_dataset(self):
		df = self.df
		folder = self.folder
		df.to_csv(folder + "TCGA_Breast_ALL_Proteome.itraq_NORM.tsv",sep = "\t")

	def translate_data_breast(self):
		data_breast = self.data_breast
		protein_list = list(data_breast.index)
		trans_df = mg.querymany(protein_list,
								scopes = 'entrezgene,ensembl.gene,symbol,reporter,accession,uniprot,name,alias,summary',
								fields = 'symbol,entrezgene,uniprot,summary,alias,go,pathway,name,homologene,homologene.id,ensembl.gene,ensembl.protein',
								species = "human",as_dataframe = True,returnall = True)
		trans_dict = trans_df["out"].to_dict()
		return trans_dict

	def map_information_breast(self):
		data_breast = self.data_breast
		trans_dict = self.trans_dict
		folder = self.folder

		protein_list = list(data_breast.index)
		entrezList = list()
		ensembl_gene_list = list()
		ensembl_protein_list = list()
		nameList = list()
		aliasList = list()
		pathway_idList = list()
		pathway_nameList = list()
		uniprotList = list()
		for protein in protein_list:
			entrezList.append(trans_dict["entrezgene"][protein])
			nameList.append(trans_dict["name"][protein])
			if type(trans_dict["alias"][protein]) == list:
				aliasList.append(",".join(trans_dict["alias"][protein]))
			else:
				aliasList.append(trans_dict["alias"][protein])
			if str(trans_dict["uniprot"][protein])!= "nan":
				get_swissprot = trans_dict["uniprot"][protein].get("Swiss-Prot",0)
				get_trembl = trans_dict["uniprot"][protein].get("TrEMBL",0)
				uniList = list()
				if get_swissprot!= 0:
					if type(trans_dict["uniprot"][protein]["Swiss-Prot"]) == list:
						for item in trans_dict["uniprot"][protein]["Swiss-Prot"]:
							uniList.append(item)
					else:
						uniList.append(trans_dict["uniprot"][protein]["Swiss-Prot"])
				if get_trembl!= 0:
					if type(trans_dict["uniprot"][protein]["TrEMBL"]) == list:
						for item in trans_dict["uniprot"][protein]["TrEMBL"]:
							uniList.append(item)
					else:
						uniList.append(trans_dict["uniprot"][protein]["TrEMBL"])
				uniprotList.append(",".join(uniList))
			else:
				uniprotList.append("")
			if str(trans_dict["pathway"][protein])!= "nan":
				get_pathway = trans_dict["pathway"][protein].get("reactome",0)
				id_patList = list()
				name_patList = list()
				if get_pathway!= 0:
					if type(trans_dict["pathway"][protein]["reactome"]) == list:
						for item in trans_dict["pathway"][protein]["reactome"]:
							id_patList.append(item["id"])
							name_patList.append(item["name"])
					else:
						id_patList.append(trans_dict["pathway"][protein]["reactome"]["id"])
						name_patList.append(trans_dict["pathway"][protein]["reactome"]["name"])
				pathway_idList.append("//".join(id_patList))
				pathway_nameList.append("//".join(name_patList))
			else:
				pathway_idList.append("")
				pathway_nameList.append("")
			if str(trans_dict["ensembl"][protein])!= "nan":
				e_genes = list()
				e_proteins = list()
				try:
					get_gene = trans_dict["ensembl"][protein].get("gene",0)
					get_protein = trans_dict["ensembl"][protein].get("protein",0)
					if get_gene!= 0:
						if type(get_gene) == list:
							for item in get_gene:
								e_genes.append(item)
						else:
							e_genes.append(get_gene)
					if get_protein!= 0:
						if type(get_protein) == list:
							for item in get_protein:
								e_proteins.append(item)
						else:
							e_proteins.append(get_protein)
				except:
					get_gene = trans_dict["ensembl"][protein][0].get("gene",0)
					get_protein = trans_dict["ensembl"][protein][0].get("protein",0)
					if get_gene!= 0:
						if type(get_gene) == list:
							for item in get_gene:
								e_genes.append(item)
						else:
							e_genes.append(get_gene)
					if get_protein!= 0:
						if type(get_protein) == list:
							for item in get_protein:
								e_proteins.append(item)
						else:
							e_proteins.append(get_protein)
				ensembl_gene_list.append(",".join(e_genes))
				ensembl_protein_list.append(",".join(e_proteins))	
			else:
				ensembl_gene_list.append("")		
				ensembl_protein_list.append("")

		data_breast["entrez"] = pd.Series(entrezList,index = data_breast.index)
		data_breast["name"] = pd.Series(nameList,index = data_breast.index)
		data_breast["alias"] = pd.Series(aliasList,index = data_breast.index)
		data_breast["uniprot"] = pd.Series(uniprotList,index = data_breast.index)
		data_breast["pathway_id"] = pd.Series(pathway_idList,index = data_breast.index)
		data_breast["pathway_name"] = pd.Series(pathway_nameList,index = data_breast.index)
		data_breast["ensembl.gene"] = pd.Series(ensembl_gene_list,index = data_breast.index)
		data_breast["ensembl.protein"] = pd.Series(ensembl_protein_list,index = data_breast.index)

		proteinList = list(data_breast.index)
		complex_idList = list()
		complex_nameList = list()
		for protein in proteinList:
			complex_ids = list()
			complex_names = list()
			for complexID in complexDict:
				humanSubunits = complexDict[complexID]["humanGeneNames"]
				if protein in humanSubunits:
					complex_ids.append(complexID)
					complex_names.append(complexDict[complexID]["altName"][0])
			complex_idList.append(",".join(complex_ids))
			complex_nameList.append(",".join(complex_names))
		data_breast["complex_formal_ID"] = pd.Series(complex_idList,index = data_breast.index)
		data_breast["complexName"] = pd.Series(complex_nameList,index = data_breast.index)
		fileName = "TCGA_Breast_ALL_Proteome.itraq_NORM_INFO.tsv"
		data_breast.to_csv(folder + fileName,sep = "\t")
		return data_breast

class step1_wider_mapping(object):
	def __init__(self,fileData,proteinID,**kwargs):
		self.fileFolder = kwargs.get('folder','PATH')
		self.species = kwargs.get("species","human")
		self.outputFileName = kwargs.get("outputFileName","")

		self.fileData = fileData
		self.fileData.index = self.fileData[proteinID]
		self.fileData = self.fileData.drop([proteinID],1)

		print("identify_index_type")
		self.trans_dict = step1_wider_mapping.identify_index_type(fileFolder, data)

		print("map_information_onto_fileData")
		if self.species == "human":
			self.mappedData = self.human_map_information_onto_fileData()
		if self.species == "mouse":
			self.mappedData = self.mouse_map_information_onto_fileData()

	@staticmethod
	def identify_index_type(folder, data):

		trans_df = mg.querymany(list(data.index),
								scopes = 'ensembl.gene,ensembl.protein,symbol,reporter,accession',
								fields = 'all',
								species = "human,mouse",
								as_dataframe = True,
								returnall = True)

		trans_df = trans_df["out"]
		trans_dict = trans_df.to_dict()
		return trans_dict

	@staticmethod
	def join_to_txt(lst):
		return ''.join(lst)

	def get_pathways(self,uni,category,pathwayList):
		trans_dict = self.trans_dict

		if str(trans_dict["pathway"][uni])!= "nan":
			getPathway = trans_dict["pathway"][uni].get(category,0)
			if getPathway!= 0:
				try:
					keys = trans_dict["pathway"][uni][category].keys()
					t = [trans_dict["pathway"][uni][category]]
				except:
					t = trans_dict["pathway"][uni][category]
				tempList = list()
				for item in t:
					try:
						tempList.append(step1_wider_mapping.join_to_txt(['id:', str(item['id']), ',name:', str(item['name']).split('\\')[0]]))
					except:
						tempList.append(step1_wider_mapping.join_to_txt(['id:', str(item['id']), ',name:', str(item["name"][0:30])]))
				pathwayList.append("/".join(tempList))
			else:
				pathwayList.append("")	
		else:
			pathwayList.append("")
		return pathwayList

	def get_refSeqs(self,uni,refSeq_geneList,refSeq_proteinList,refSeq_rnaList,refSeq_translationList):
		trans_dict = self.trans_dict

		if str(trans_dict["refseq"][uni])!= "nan":
			if trans_dict["refseq"][uni].get("genomic",0)!= 0:
				t = trans_dict["refseq"][uni]["genomic"]
				if type(t)!= list:
					refSeq_geneList.append(t)
				else:
					refSeq_geneList.append(",".join(t))
			else:
				refSeq_geneList.append("")
			if trans_dict["refseq"][uni].get("protein",0)!= 0:
				t = trans_dict["refseq"][uni]["protein"]
				if type(t)!= list:
					refSeq_proteinList.append(t)
				else:
					refSeq_proteinList.append(",".join(t))
			else:
				refSeq_proteinList.append("")
			if trans_dict["refseq"][uni].get("rna",0)!= 0:
				t = trans_dict["refseq"][uni]["rna"]
				if type(t)!= list:
					refSeq_rnaList.append(t)
				else:
					refSeq_rnaList.append(",".join(t))
			else:
				refSeq_rnaList.append("")
			if trans_dict["refseq"][uni].get("translation",0)!= 0:
				try:
					keys = trans_dict["refseq"][uni]["translation"].keys()
					t = [trans_dict["refseq"][uni]["translation"]]
				except:
					t = trans_dict["refseq"][uni]["translation"]
				tempList = list()
				for item in t:
					tempList.append(step1_wider_mapping.join_to_txt(['protein:', str(item['protein']), ',rna:', str(item['rna'])]))
				refSeq_translationList.append(",".join(tempList))
			else:
				refSeq_translationList.append("")
		else:
			refSeq_geneList.append("")
			refSeq_proteinList.append("")
			refSeq_rnaList.append("")
			refSeq_translationList.append("")
		return refSeq_geneList,refSeq_proteinList,refSeq_rnaList,refSeq_translationList

	def get_ensembl_information(self,uni):
		trans_dict = self.trans_dict

		geneList = list()
		proteinList = list()
		transcriptList = list()
		translationList = list()
		if str(trans_dict["ensembl"][uni])!= "nan":
			try:
				keys = trans_dict["ensembl"][uni].keys()
				t = [trans_dict["ensembl"][uni]]
			except:
				t = trans_dict["ensembl"][uni]
			for item in t:
				if item.get("gene",0)!= 0:
					if type(item["gene"])!= list:
						geneList.append(item["gene"])
					else:
						geneList.append(",".join(item["gene"]))
				if item.get("protein",0)!= 0:
					if type(item["protein"])!= list:
						proteinList.append(item["protein"])
					else:
						proteinList.append(",".join(item["protein"]))
				if item.get("transcript",0)!= 0:
					if type(item["transcript"])!= list:
						transcriptList.append(item["transcript"])
					else:
						transcriptList.append(",".join(item["transcript"]))
				if item.get("translation",0)!= 0:
					translation = ""
					for trans in item["translation"]:
						translation+= trans["protein"] + ":" + trans["rna"] + ","
					translationList.append(translation[:-1])
		else:
			geneList.append("")
			proteinList.append("")
			transcriptList.append("")
			translationList.append("")
		return geneList,proteinList,transcriptList,translationList

	def get_genomic_position(self,uni):
		trans_dict = self.trans_dict

		if self.species == "human":
			cat = "hg19"
		elif self.species == "mouse":
			cat = "mm9"

		genomicPosList = list()
		genomicPos_hg19_List = list()
		if str(trans_dict["genomic_pos"][uni])!= "nan" and str(trans_dict["genomic_pos_" + cat][uni])!= "nan":
			try:
				keys = trans_dict["genomic_pos"][uni].keys()
				t = [trans_dict["genomic_pos"][uni]]
			except:
				t = trans_dict["genomic_pos"][uni]
			if self.species!= "fly":
				try:
					keys = trans_dict["genomic_pos_" + cat][uni].keys()
					t_hg19 = [trans_dict["genomic_pos_" + cat][uni]]
				except:
					t_hg19 = trans_dict["genomic_pos_" + cat][uni]
			for item in t:
				temp = step1_wider_mapping.join_to_txt(['chr:', str(item['chr']), ',start:', str(item['start']), ',end:', str(item['end']), ',strand:', str(item['strand'])])
				genomicPosList.append(temp)
			if self.species!= "fly":
				try:
					for item in t_hg19:
						temp = step1_wider_mapping.join_to_txt(['chr:', str(item['chr']), ',start:', str(item['start']), ',end:', str(item['end']), ',strand:', str(item['strand'])])
						genomicPos_hg19_List.append(temp)
				except:
					raise Exception(t_hg19)
		else:
			genomicPosList.append("")
			genomicPos_hg19_List.append("")
		return genomicPosList,genomicPos_hg19_List

	def get_interpro_domains(self,uni):
		trans_dict = self.trans_dict

		if str(trans_dict["interpro"][uni])!= "nan":
			try:
				keys = trans_dict["interpro"][uni].keys()
				t = [trans_dict["interpro"][uni]]
			except:
				t = trans_dict["interpro"][uni]
			tempList = list()
			for item in t:
				tempList.append(step1_wider_mapping.join_to_txt(['desc:', item['desc'], ',id:', item['id'], ',short_desc:', item['short_desc']]))
			temp = "/".join(tempList)
		else:
			temp = ""
		return temp

	def get_GO(self,uni):
		trans_dict = self.trans_dict
		if str(trans_dict["go"][uni])! = "nan":
			get_bp = trans_dict["go"][uni].get("BP",0)
			if get_bp!= 0:
				tempList = list()
				try:
					keys = trans_dict["go"][uni]["BP"].keys()
					transList = [trans_dict["go"][uni]["BP"]]
				except:
					transList = trans_dict["go"][uni]["BP"]
				for t in transList:
					tempList.append(step1_wider_mapping.join_to_txt(['evidence:', t['evidence'], ',id:', t['id'], ',term:', t['term']]))
				bp = "/".join(tempList)
			else:
				bp = ""
			get_mf = trans_dict["go"][uni].get("MF",0)
			if get_mf!= 0:	
				tempList = list()
				try:
					keys = trans_dict["go"][uni]["MF"].keys()
					transList = [trans_dict["go"][uni]["MF"]]
				except:
					transList = trans_dict["go"][uni]["MF"]
				for t in transList:
					tempList.append(step1_wider_mapping.join_to_txt(['evidence:', t['evidence'], ',id:', t['id'], ',term:', t['term']]))
				mf = "/".join(tempList)
			else:
				mf = ""
			get_cc = trans_dict["go"][uni].get("CC",0)
			if get_cc!= 0:		
				tempList = list()
				try:
					keys = trans_dict["go"][uni]["CC"].keys()
					transList = [trans_dict["go"][uni]["CC"]]
				except:
					transList = trans_dict["go"][uni]["CC"]
				for t in transList:
					tempList.append(step1_wider_mapping.join_to_txt(['evidence:', t['evidence'], ',id:', t['id'], ',term:', t['term']]))
				cc = "/".join(tempList)
			else:
				cc = ""
		else:
			bp = "";cc = "";mf = ""
		return bp,mf,cc

	def get_exons(self,uni):
		if self.species == "human":
			cat = "hg19"
		if self.species == "mouse":
			cat = "mm9"
		trans_dict = self.trans_dict
		if str(trans_dict["exons"][uni])! = "nan":
			try:
				keys = trans_dict["exons"][uni].keys()
				ex = [trans_dict["exons"][uni]]
			except:
				ex = trans_dict["exons"][uni]
			tempList = list()
			for item in ex:
				txt1 = step1_wider_mapping.join_to_txt(['transcript:', item['transcript'],
														',chr:', str(item['chr']), ',cdsstart:',
														str(item['cdsstart']), ',cdsend:', str(item['cdsend'])])
				txt1 = txt1 + ','
				txt2 = step1_wider_mapping.join_to_txt(['strand:', str(item['strand']), ',txstart:',
														str(item['txstart']), ',txend:', str(item['txend']),
														',positions:'])
				posList = list()
				for pos in item["position"]:
					posList.append(str(pos[0]) + "-" + str(pos[1]))
				txt = txt1 + txt2 + ";".join(posList)
				tempList.append(txt)
			exons = "/".join(tempList)
		else:
			exons = ""
		exons_hg19 = ""

		if self.species!= "fly":
			exons_hg19 = ""
			if str(trans_dict["exons_" + cat][uni])!= "nan":
				try:
					keys = trans_dict["exons_" + cat][uni].keys()
					ex = [trans_dict["exons_" + cat][uni]]
				except:
					ex = trans_dict["exons"][uni]
				ex = trans_dict["exons_" + cat][uni]
				tempList = list()
				for item in ex:
					txt1 = step1_wider_mapping.join_to_txt(['transcript:', item['transcript'],
															',chr:', str(item['chr']), ',cdsstart:',
															str(item['cdsstart']), ',cdsend:', str(item['cdsend'])])
					txt1 = txt1 + ','
					txt2 = step1_wider_mapping.join_to_txt(['strand:', str(item['strand']), ',txstart:',
															str(item['txstart']), ',txend:', str(item['txend']),
															',positions:'])
					posList = list()
					for pos in item["position"]:
						posList.append(str(pos[0]) + "-" + str(pos[1]))
					txt = txt1 + txt2 + ";".join(posList)
					tempList.append(txt)
				exons_hg19 = "/".join(tempList)
			else:
				exons_hg19 = ""
		return exons,exons_hg19

	def get_uniprot(self,uni):
		trans_dict = self.trans_dict
		if str(trans_dict["uniprot"][uni])! = "nan":
			tempList = list()
			get_sw = trans_dict["uniprot"][uni].get("Swiss-Prot",0)
			if get_sw!= 0:
				sw_uniList = trans_dict["uniprot"][uni]["Swiss-Prot"]
			get_tr = trans_dict["uniprot"][uni].get("TrEMBL",0)
			if get_tr!= 0:
				tr_uniList = trans_dict["uniprot"][uni]["TrEMBL"]
			if get_sw!= 0:
				if len(sw_uniList)>1 and len(sw_uniList[0])>= 6:
					for item in sw_uniList:
						tempList.append(item)
				else:
					tempList.append(sw_uniList)
			if get_tr!= 0:
				if len(tr_uniList)>1 and len(tr_uniList[0])>= 6:
					for item in tr_uniList:
						tempList.append(item)
				else:
					tempList.append(tr_uniList)	
			uniprot = ",".join(tempList)
		else:
			uniprot = ""
		return uniprot

	def get_pubmed(self,uni):
		trans_dict = self.trans_dict
		if str(trans_dict["generif"][uni])!= "nan":
			try:
				keys = trans_dict["generif"][uni].keys()
				transList = [trans_dict["generif"][uni]]
			except:
				transList = trans_dict["generif"][uni]
			tempList = list()
			for t in transList:
				tempList.append(str(t["pubmed"]))#+":"+t["text"])	
			pubmed = ",".join(tempList)
		else:
			pubmed = ""
		return pubmed

	def get_homologenes(self,uni):
		trans_dict = self.trans_dict
		t = trans_dict["homologene"][uni]
		if str(t)!= "nan":
			temp = "homologene:" + str(t["id"]) + ","
			for item in t["genes"]:
				temp+ = "taxID:" + str(item[0]) + ":" + str(item[1]) + "(entrez),"
			temp = temp[:-1]
		else:
			temp = ""
		return temp

	def get_more_info(self,uni):
		trans_dict = self.trans_dict
		if str(trans_dict["alias"][uni])!= "nan":
			if type(trans_dict["alias"][uni])!= list:
				transList = [trans_dict["alias"][uni]]
			else:
				transList = trans_dict["alias"][uni]
			alias = ",".join(transList)
		else:
			alias = ""
		if self.species!= "fly":
			if str(trans_dict["ipi"][uni])!= "nan":
				if type(trans_dict["ipi"][uni])!= list:
					transList = [trans_dict["ipi"][uni]]
				else:
					transList = trans_dict["ipi"][uni]
				ipi = ",".join(transList)
			else:
				ipi = ""
		else:
			ipi = ""
		if str(trans_dict["other_names"][uni])!= "nan":
			if type(trans_dict["other_names"][uni])!= list:
				transList = [trans_dict["other_names"][uni]]
			else:
				transList = trans_dict["other_names"][uni]
			other_names = ",".join(transList)
		else:
			other_names = ""
		if str(trans_dict["pdb"][uni])!= "nan":
			if type(trans_dict["pdb"][uni])!= list:
				transList = [trans_dict["pdb"][uni]]
			else:
				transList = trans_dict["pdb"][uni]
			pdb = ",".join(transList)
		else:
			pdb = ""
		if str(trans_dict["pfam"][uni])!= "nan":
			if type(trans_dict["pfam"][uni])!= list:
				transList = [trans_dict["pfam"][uni]]
			else:
				transList = trans_dict["pfam"][uni]
			pfam = ",".join(transList)
		else:
			pfam = ""
		if str(trans_dict["prosite"][uni])!= "nan":
			if type(trans_dict["prosite"][uni])!= list:
				transList = [trans_dict["prosite"][uni]]
			else:
				transList = trans_dict["prosite"][uni]
			prosite = ",".join(transList)
		else:
			prosite = ""
		if str(trans_dict["unigene"][uni])!= "nan":
			if type(trans_dict["unigene"][uni])!= list:
				transList = [trans_dict["unigene"][uni]]
			else:
				transList = trans_dict["unigene"][uni]
			unigene = ",".join(transList)
		else:
			unigene = ""
		return {"alias":alias,
				"ipi":ipi,
				"other_names":other_names,
				"pdb":pdb,
				"pfam":pfam,
				"prosite":prosite,
				"unigene":unigene}

	def human_map_information_onto_fileData(self):
		fileData = self.fileData
		trans_dict = self.trans_dict
		fileFolder = self.fileFolder

		entrez = list()
		hgnc_list = list()
		mim_list = list()
		vega_list = list()
		aliasList = list()
		ec_list = list()
		pubmedList = list()
		geneTypeList = list()
		ipiList = list()
		nameList = list()
		homologeneList = list()
		mapLocationList = list()
		otherNameList = list()
		exons = list();exons_hg19 = list()
		ensembl_protein = list();ensembl_gene = list();ensembl_transcript = list();ensembl_translation = list()
		genomic_pos = list();genomic_pos_hg19 = list()
		go_bpList = list();go_mfList = list();go_ccList = list()
		pathway_biocarta_List = list();pathway_humancyc_List = list();pathway_kegg_List = list()
		pathway_pharmgkb_List = list();pathway_reactome_List = list();pathway_smpdb_List = list();pathway_wiki_List = list()
		refSeq_proteinList = list();refSeq_geneList = list();refSeq_rnaList = list();refSeq_translationList = list()
		pdbList = list();pfamList = list();prositeList = list();interproList = list();pharmgkbList = list()
		summaryList = list();symbolList = list()
		unigeneList = list();uniprotList = list()

		for u,uni in enumerate(list(fileData.index)):
			print(u)
			getKey = trans_dict["uniprot"].get(uni,0)
			if getKey!= 0:
				hgnc_list.append(trans_dict["HGNC"][uni])
				mim_list.append(trans_dict["MIM"][uni])
				vega_list.append(trans_dict["Vega"][uni])
				nameList.append("/".join(str(trans_dict["name"][uni]).split(",")))
				ec_list.append(trans_dict["ec"][uni])
				mapLocationList.append(",".join(str(trans_dict["map_location"][uni]).split(",")))
				pharmgkbList.append(",".join(str(trans_dict["pharmgkb"][uni]).split(",")))

				geneList,proteinList,transcriptList,translationList = self.get_ensembl_information(uni)
				ensembl_gene.append("/".join(geneList))
				ensembl_protein.append("/".join(proteinList))
				ensembl_transcript.append("/".join(transcriptList))
				ensembl_translation.append("/".join(translationList))

				entrez.append(trans_dict["entrezgene"][uni])
				pubmed = self.get_pubmed(uni)
				pubmedList.append(pubmed)

				posList,pos_hg19_List = self.get_genomic_position(uni)
				genomic_pos.append("/".join(posList))
				genomic_pos_hg19.append("/".join(pos_hg19_List))

				bp,mf,cc = self.get_GO(uni)
				go_bpList.append(bp)
				go_mfList.append(mf)
				go_ccList.append(cc)

				homologene = self.get_homologenes(uni)
				homologeneList.append(homologene)

				infoDict = self.get_more_info(uni)
				aliasList.append(infoDict["alias"])
				ipiList.append(infoDict["ipi"])
				otherNameList.append(infoDict["other_names"])
				pdbList.append(infoDict["pdb"])
				pfamList.append(infoDict["pfam"])
				prositeList.append(infoDict["prosite"])
				unigeneList.append(infoDict["unigene"])

				pathway_biocarta_List = self.get_pathways(uni,"biocarta",pathway_biocarta_List)
				pathway_humancyc_List = self.get_pathways(uni,"humancyc",pathway_humancyc_List)
				pathway_kegg_List = self.get_pathways(uni,"kegg",pathway_kegg_List)
				pathway_pharmgkb_List = self.get_pathways(uni,"pharmgkb",pathway_pharmgkb_List)
				pathway_reactome_List = self.get_pathways(uni,"reactome",pathway_reactome_List)
				pathway_smpdb_List = self.get_pathways(uni,"smpdb",pathway_smpdb_List)
				pathway_wiki_List = self.get_pathways(uni,"wikipathways",pathway_wiki_List)

				refSeq_geneList,refSeq_proteinList,refSeq_rnaList,refSeq_translationList = self.get_refSeqs(uni,refSeq_geneList,refSeq_proteinList,refSeq_rnaList,refSeq_translationList)
				interpro = self.get_interpro_domains(uni)
				interproList.append(interpro)

				summaryList.append(trans_dict["summary"][uni])
				symbolList.append(trans_dict["symbol"][uni])
				geneTypeList.append(trans_dict["type_of_gene"][uni])
				
				ex,ex_hg19 = self.get_exons(uni)
				exons.append(ex)
				exons_hg19.append(ex_hg19)

				uniprot = self.get_uniprot(uni)
				uniprotList.append(uniprot)

		fileData["HGNC"] = pd.Series(hgnc_list,index = fileData.index)
		fileData["MIM"] = pd.Series(mim_list,index = fileData.index)
		fileData["VEGA"] = pd.Series(vega_list,index = fileData.index)
		fileData["ALIAS"] = pd.Series(aliasList,index = fileData.index)
		fileData["EC"] = pd.Series(ec_list,index = fileData.index)
		fileData["ENSEMBL_PROTEIN"] = pd.Series(ensembl_protein,index = fileData.index)
		fileData["ENSEMBL_GENE"] = pd.Series(ensembl_gene,index = fileData.index)
		fileData["ENSEMBL_TRANSCRIPT"] = pd.Series(ensembl_transcript,index = fileData.index)
		fileData["ENSEMBL_TRANSLATION"] = pd.Series(ensembl_translation,index = fileData.index)
		fileData["ENTREZ"] = pd.Series(entrez,index = fileData.index)
		fileData["EXONS"] = pd.Series(exons,index = fileData.index)
		fileData["EXONS_HG19"] = pd.Series(exons_hg19,index = fileData.index)
		fileData["PUBMED"] = pd.Series(pubmedList,index = fileData.index)
		fileData["GENOMIC_POS"] = pd.Series(genomic_pos,index = fileData.index)
		fileData["GENOMIC_POS_HG19"] = pd.Series(genomic_pos_hg19,index = fileData.index)
		fileData["GO_BP"] = pd.Series(go_bpList,index = fileData.index)
		fileData["GO_MF"] = pd.Series(go_mfList,index = fileData.index)
		fileData["GO_CC"] = pd.Series(go_ccList,index = fileData.index)
		fileData["HOMOLOGENE"] = pd.Series(homologeneList,index = fileData.index)
		fileData["INTERPRO"] = pd.Series(interproList,index = fileData.index)
		fileData["IPI"] = pd.Series(ipiList,index = fileData.index)
		fileData["NAME"] = pd.Series(nameList,index = fileData.index)
		fileData["MAP_LOCATION"] = pd.Series(mapLocationList,index = fileData.index)
		fileData["OTHER_NAME"] = pd.Series(otherNameList,index = fileData.index)
		fileData["PATHWAY_BIOCARTA"] = pd.Series(pathway_biocarta_List,index = fileData.index)
		fileData["PATHWAY_HUMANCYC"] = pd.Series(pathway_humancyc_List,index = fileData.index)
		fileData["PATHWAY_KEGG"] = pd.Series(pathway_kegg_List,index = fileData.index)
		fileData["PATHWAY_PHARMGKB"] = pd.Series(pathway_pharmgkb_List,index = fileData.index)
		fileData["PATHWAY_REACTOME"] = pd.Series(pathway_reactome_List,index = fileData.index)
		fileData["PATHWAY_SMPDB"] = pd.Series(pathway_smpdb_List,index = fileData.index)
		fileData["PATHWAY_WIKI"] = pd.Series(pathway_wiki_List,index = fileData.index)
		fileData["PDB"] = pd.Series(pdbList,index = fileData.index)
		fileData["PFAM"] = pd.Series(pfamList,index = fileData.index)
		fileData["PHARMGKB"] = pd.Series(pharmgkbList,index = fileData.index)
		fileData["PROSITE"] = pd.Series(prositeList,index = fileData.index)
		fileData["REFSEQ_PROTEIN"] = pd.Series(refSeq_proteinList,index = fileData.index)
		fileData["REFSEQ_GENE"] = pd.Series(refSeq_geneList,index = fileData.index)
		fileData["REFSEQ_RNA"] = pd.Series(refSeq_rnaList,index = fileData.index)
		fileData["REFSEQ_TRANSLATION"] = pd.Series(refSeq_translationList,index = fileData.index)
		fileData["SUMMARY"] = pd.Series(summaryList,index = fileData.index)
		fileData["SYMBOL"] = pd.Series(symbolList,index = fileData.index)
		fileData["GENE_TYPE"] = pd.Series(geneTypeList,index = fileData.index)
		fileData["UNIGENE"] = pd.Series(unigeneList,index = fileData.index)
		fileData["UNIPROT"] = pd.Series(uniprotList,index = fileData.index)

		fileData.to_csv(fileFolder + time.strftime("%Y%m%d") + "_" + self.outputFileName + "_ALL_INFORMATION.txt",
						sep = "\t")
		return fileData

	def mouse_map_information_onto_fileData(self):
		fileData = self.fileData
		trans_dict = self.trans_dict
		fileFolder = self.fileFolder

		entrez = list()
		hgnc_list = list()
		mim_list = list()
		vega_list = list()
		aliasList = list()
		ec_list = list()
		pubmedList = list()
		geneTypeList = list()
		ipiList = list()
		nameList = list()
		homologeneList = list()
		mapLocationList = list()
		otherNameList = list()
		exons = list();exons_hg19 = list()
		ensembl_protein = list();ensembl_gene = list();ensembl_transcript = list();ensembl_translation = list()
		genomic_pos = list();genomic_pos_hg19 = list()
		go_bpList = list();go_mfList = list();go_ccList = list()
		pathway_biocarta_List = list();pathway_humancyc_List = list();pathway_kegg_List = list()
		pathway_pharmgkb_List = list();pathway_reactome_List = list();pathway_smpdb_List = list();pathway_wiki_List = list()
		refSeq_proteinList = list();refSeq_geneList = list();refSeq_rnaList = list();refSeq_translationList = list()
		pdbList = list();pfamList = list();prositeList = list();interproList = list();pharmgkbList = list()
		summaryList = list();symbolList = list()
		unigeneList = list();uniprotList = list()

		for u,uni in enumerate(list(fileData.index)):
			print(u)
			getKey = trans_dict["uniprot"].get(uni,0)
			if getKey!= 0:
				mim_list.append(trans_dict["MGI"][uni])#MIM
				vega_list.append(trans_dict["Vega"][uni])
				nameList.append("/".join(str(trans_dict["name"][uni]).split(",")))
				ec_list.append(trans_dict["ec"][uni])
				mapLocationList.append(",".join(str(trans_dict["map_location"][uni]).split(",")))

				geneList,proteinList,transcriptList,translationList = self.get_ensembl_information(uni)
				ensembl_gene.append("/".join(geneList))
				ensembl_protein.append("/".join(proteinList))
				ensembl_transcript.append("/".join(transcriptList))
				ensembl_translation.append("/".join(translationList))

				entrez.append(trans_dict["entrezgene"][uni])
				pubmed = self.get_pubmed(uni)
				pubmedList.append(pubmed)

				posList,pos_hg19_List = self.get_genomic_position(uni)
				genomic_pos.append("/".join(posList))
				genomic_pos_hg19.append("/".join(pos_hg19_List))

				bp,mf,cc = self.get_GO(uni)
				go_bpList.append(bp)
				go_mfList.append(mf)
				go_ccList.append(cc)

				homologene = self.get_homologenes(uni)
				homologeneList.append(homologene)

				infoDict = self.get_more_info(uni)
				aliasList.append(infoDict["alias"])
				ipiList.append(infoDict["ipi"])
				otherNameList.append(infoDict["other_names"])
				pdbList.append(infoDict["pdb"])
				pfamList.append(infoDict["pfam"])
				prositeList.append(infoDict["prosite"])
				unigeneList.append(infoDict["unigene"])

				pathway_biocarta_List = self.get_pathways(uni,"biocarta",pathway_biocarta_List)
				pathway_humancyc_List = self.get_pathways(uni,"humancyc",pathway_humancyc_List)
				pathway_kegg_List = self.get_pathways(uni,"kegg",pathway_kegg_List)
				pathway_pharmgkb_List = self.get_pathways(uni,"pharmgkb",pathway_pharmgkb_List)
				pathway_reactome_List = self.get_pathways(uni,"reactome",pathway_reactome_List)
				pathway_smpdb_List = self.get_pathways(uni,"smpdb",pathway_smpdb_List)
				pathway_wiki_List = self.get_pathways(uni,"wikipathways",pathway_wiki_List)

				refSeq_geneList,refSeq_proteinList,refSeq_rnaList,refSeq_translationList = self.get_refSeqs(uni,refSeq_geneList,refSeq_proteinList,refSeq_rnaList,refSeq_translationList)
				interpro = self.get_interpro_domains(uni)
				interproList.append(interpro)

				summaryList.append(trans_dict["summary"][uni])
				symbolList.append(trans_dict["symbol"][uni])
				geneTypeList.append(trans_dict["type_of_gene"][uni])
				
				ex,ex_hg19 = self.get_exons(uni)
				exons.append(ex)
				exons_hg19.append(ex_hg19)

				uniprot = self.get_uniprot(uni)
				uniprotList.append(uniprot)

		fileData["MGI"] = pd.Series(mim_list,index = fileData.index)
		fileData["VEGA"] = pd.Series(vega_list,index = fileData.index)
		fileData["ALIAS"] = pd.Series(aliasList,index = fileData.index)
		fileData["EC"] = pd.Series(ec_list,index = fileData.index)
		fileData["ENSEMBL_PROTEIN"] = pd.Series(ensembl_protein,index = fileData.index)
		fileData["ENSEMBL_GENE"] = pd.Series(ensembl_gene,index = fileData.index)
		fileData["ENSEMBL_TRANSCRIPT"] = pd.Series(ensembl_transcript,index = fileData.index)
		fileData["ENSEMBL_TRANSLATION"] = pd.Series(ensembl_translation,index = fileData.index)
		fileData["ENTREZ"] = pd.Series(entrez,index = fileData.index)
		fileData["EXONS"] = pd.Series(exons,index = fileData.index)
		fileData["EXONS_MM9"] = pd.Series(exons_hg19,index = fileData.index)#HG19
		fileData["PUBMED"] = pd.Series(pubmedList,index = fileData.index)
		fileData["GENOMIC_POS"] = pd.Series(genomic_pos,index = fileData.index)
		fileData["GENOMIC_POS_MM9"] = pd.Series(genomic_pos_hg19,index = fileData.index)#HG19
		fileData["GO_BP"] = pd.Series(go_bpList,index = fileData.index)
		fileData["GO_MF"] = pd.Series(go_mfList,index = fileData.index)
		fileData["GO_CC"] = pd.Series(go_ccList,index = fileData.index)
		fileData["HOMOLOGENE"] = pd.Series(homologeneList,index = fileData.index)
		fileData["INTERPRO"] = pd.Series(interproList,index = fileData.index)
		fileData["IPI"] = pd.Series(ipiList,index = fileData.index)
		fileData["NAME"] = pd.Series(nameList,index = fileData.index)
		fileData["MAP_LOCATION"] = pd.Series(mapLocationList,index = fileData.index)
		fileData["OTHER_NAME"] = pd.Series(otherNameList,index = fileData.index)
		fileData["PATHWAY_BIOCARTA"] = pd.Series(pathway_biocarta_List,index = fileData.index)
		fileData["PATHWAY_HUMANCYC"] = pd.Series(pathway_humancyc_List,index = fileData.index)
		fileData["PATHWAY_KEGG"] = pd.Series(pathway_kegg_List,index = fileData.index)
		fileData["PATHWAY_PHARMGKB"] = pd.Series(pathway_pharmgkb_List,index = fileData.index)
		fileData["PATHWAY_REACTOME"] = pd.Series(pathway_reactome_List,index = fileData.index)
		fileData["PATHWAY_SMPDB"] = pd.Series(pathway_smpdb_List,index = fileData.index)
		fileData["PATHWAY_WIKI"] = pd.Series(pathway_wiki_List,index = fileData.index)
		fileData["PDB"] = pd.Series(pdbList,index = fileData.index)
		fileData["PFAM"] = pd.Series(pfamList,index = fileData.index)
		fileData["PROSITE"] = pd.Series(prositeList,index = fileData.index)
		fileData["REFSEQ_PROTEIN"] = pd.Series(refSeq_proteinList,index = fileData.index)
		fileData["REFSEQ_GENE"] = pd.Series(refSeq_geneList,index = fileData.index)
		fileData["REFSEQ_RNA"] = pd.Series(refSeq_rnaList,index = fileData.index)
		fileData["REFSEQ_TRANSLATION"] = pd.Series(refSeq_translationList,index = fileData.index)
		fileData["SUMMARY"] = pd.Series(summaryList,index = fileData.index)
		fileData["SYMBOL"] = pd.Series(symbolList,index = fileData.index)
		fileData["GENE_TYPE"] = pd.Series(geneTypeList,index = fileData.index)
		fileData["UNIGENE"] = pd.Series(unigeneList,index = fileData.index)
		fileData["UNIPROT"] = pd.Series(uniprotList,index = fileData.index)

		fileData.to_csv(fileFolder + time.strftime("%Y%m%d") + "_" + self.outputFileName + "_ALL_INFORMATION.txt", 
						sep = "\t")
		return fileData

class step1_module_mapping:
	@staticmethod
	def execute(data, **kwargs):
		folder  =  kwargs.get('folder','PATH')
		species  =  kwargs.get('species','mouse')
		output_folder  =  kwargs.get('output_folder','PATH')
		output_fileName  =  kwargs.get('output_fileName','NAME')
		user_mappings  =  kwargs.get('user_mappings', ['complexes',
													   'reactome',
													   'compartment',
													   'essential',
													   'housekeeping',
													   'string',
													   'chromosome'])

		data  =  step1_mapping.relevant_columns_reformat(data)

		if 'reactome' in user_mappings:
			print('map_reactome_pathways')
			data  =  step1_mapping.map_reactome_pathways(data, species)

		if 'compartment' in user_mappings:
			print('map_cellular_compartment')
			data  =  step1_mapping.map_cellular_compartment(data, species, folder)

		if 'complexes' in user_mappings:
			print('map_complexes')
			data  =  step1_mapping.map_complexes(data, species)

		if 'housekeeping' in user_mappings:
			print('map_housekeeping_genes')
			data  =  step1_mapping.map_housekeeping_genes(data, folder)

		if 'essential' in user_mappings:
			print('map_essential_genes')
			data  =  step1_mapping.map_essential_genes(data, folder)

		if 'string' in user_mappings:
			print('map_string_interactions')
			data  =  step1_mapping.map_string_interactions(data, folder, species)

		if 'chromosome' in user_mappings:
			print('map_chromosome_interactions')
			data  =  step1_mapping.map_chromosome_interactions()

		print('export_dataframe')
		step1_mapping.export_dataframe(data, output_folder, output_fileName)

	@staticmethod
	def load_complex_dictionary(**kwargs):
		folder = kwargs.get('folder','PATH')

		complexDict = DataFrameAnalyzer.read_pickle(folder + 'complex_dictionary.pkl')
		return complexDict

	@staticmethod
	def relevant_columns_reformat(data):
		quant_list = ['quant_' + item for item in list(data.columns)]
		data.columns = quant_list
		return data

	@staticmethod
	def map_reactome_pathways(data, species):
		gene_list = list(set(data.index))

		m = Mapper(gene_list,
				   input = 'symbol,uniprot,entrezgene,ensembl.gene,ensembl.protein',
				   output = 'pathway.reactome',
				   species = species)

		hit_dict = m.hitDict['pathway']

		gene_list = list(data.index)
		reactome_id_list = list()
		reactome_name_list = list()
		reactome_pathway_list = list()

		for gene in gene_list:
			reactomes = []
			reactome_ids = []
			reactome_names = []
			try:
				if str(hit_dict[gene])!= 'nan':
					if type(hit_dict[gene]['reactome']) == list:
						for item in hit_dict[gene]['reactome']:
							reactomes.append(item['id'] + '::' + item['name'])
							reactome_ids.append(item['id'])
							reactome_names.append(item['name'])
					else:
						item  =  hit_dict[gene]['reactome']
						reactomes.append(item['id'] + '::' + item['name'])
						reactome_ids.append(item['id'])
						reactome_names.append(item['name'])
			except:
				pass
			reactome_id_list.append(';'.join(reactome_ids))
			reactome_name_list.append(';'.join(reactome_names))
			reactome_pathway_list.append(';'.join(reactomes))

		data['reactome_id'] = pd.Series(reactome_id_list, index = data.index)
		data['reactome_name'] = pd.Series(reactome_name_list, index = data.index)
		data['reactome'] = pd.Series(reactome_pathway_list, index = data.index)

		return data

	@staticmethod
	def get_locations(folder):
		folder = self.folder
		fileName = "subcellular_location.csv"

		locData = DataFrameAnalyzer.getFile(folder,fileName, delimiter  =  ',')
		location_list = list()
		for i,row in locData.iterrows():
			item_list = filter(lambda a:str(a)! = "nan",[row["Validated"], row["Supported"], row["Approved"]])
			if len(item_list) == 0:
				location_list.append("")
			else:
				location_list.append(";".join(item_list))
		locData["location"] = pd.Series(location_list,index  =  locData.index)
		locData.index = list(locData["Gene name"])
		locDict = locData["location"].to_dict()

		loc_set = list(locData['location'])
		lset = list()
		for l in loc_set:
			for item in l.split(';'):
				lset.append(item)
		lset = filter(lambda a:str(a)! = '',list(set(lset)))
		
		loc_name_dict = dict((e1,list()) for e1 in lset)
		for loc in lset:
			sub = locData[locData.location.str.contains(loc)]
			sub_specific = locData[locData.location =  = loc]
			loc_name_dict[loc] = {'all':list(sub.index), 
								  'specific':list(sub_specific.index)}		
		return locDict, loc_name_dict

	@staticmethod
	def map_cellular_compartment(data, species, folder):
		print('get_locations')
		locDict, loc_name_dict = step1_mapping.get_locations(folder)
		gene_list = list(set(data.index))

		m = Mapper(gene_list, 
				   input = 'symbol,uniprot,entrezgene,ensembl.gene,ensembl.protein',
				   output = 'symbol',
				   species = species)

		symbol_dict = m.hitDict['symbol']

		gene_list = list(data.index)
		location_list = list()
		for gene in gene_list:
			try:
				sym = symbol_dict[gene].upper()
				location_list.append(locDict[sym])
			except:
				location_list.append('')
		data['location_humanProteinAtlas'] = pd.Series(location_list, index = data.index)
		return data

	@staticmethod
	def map_complexes(data, species, folder):
		complexDict = step1_mapping.load_complex_dictionary(folder  =  folder)

		gene_list = list(data.index)
		complex_ids = list()
		complex_names = list()
		complex_search_list = list()
		for g,gene in enumerate(gene_list):
			gene = str(gene)
			for protein in gene.split('|'):
				com_ids, com_names = complexFacade.get_associated_complexes(protein,
									 complexDict, species = species)

			if com_ids!= '':
				complex_ids.append(com_ids)
				complex_names.append(com_names)
				com_searches = list()
				for i,n in zip(com_ids.split(';;'), com_names.split(';;')):
					com_searches.append(i + '::' + n)
				complex_search_list.append(';;'.join(com_searches))
			else:
				complex_ids.append('')
				complex_names.append('')
				complex_search_list.append('')
		data['complex_id'] = pd.Series(complex_ids, index = data.index)
		data['complex_name'] = pd.Series(complex_names, index = data.index)
		data['complex_search'] = pd.Series(complex_search_list, index = data.index)
		return data

	@staticmethod
	def map_housekeeping_genes(data, folder):
		file_name = "housekeeping_genes.txt"
		housekeeping_data = DataFrameAnalyzer.getFile(folder,file_name)
		housekeeping_genes = [item[:-1] for item in list(housekeeping_data.index)]

		gene_list = list(set(data.index))

		m = Mapper(gene_list, 
				   input = 'symbol,uniprot,entrezgene,ensembl.gene,ensembl.protein',
				   output = 'symbol',
				   species = self.species)

		symbol_dict = m.hitDict['symbol']

		gene_list = list(data.index)
		housekeeping_list = list()
		for gene in gene_list:
			try:
				sym = symbol_dict[gene]
				if sym.upper() in housekeeping_genes:
					housekeeping_list.append(sym)
				else:
					housekeeping_list.append('')
			except:
				housekeeping_list.append('')
		data['housekeeping'] = pd.Series(housekeeping_list, index = data.index)
		return data

	@staticmethod
	def map_essential_genes(data, folder):
		file_name = "essentiality_genes.txt"
		essentiality_data = DataFrameAnalyzer.getFile(folder,file_name)
		kbm7_essentiality_data = essentiality_data[essentiality_data["KBM7 adjusted p-value"]<0.05]
		k562_essentiality_data = essentiality_data[essentiality_data["K562 adjusted p-value"]<0.05]
		jiyoye_essentiality_data = essentiality_data[essentiality_data["Jiyoye adjusted p-value"]<0.05]
		raji_essentiality_data = essentiality_data[essentiality_data["Raji adjusted p-value"]<0.05]
		kbm7_essential_genes = kbm7_essentiality_data.index
		k562_essential_genes = k562_essentiality_data.index
		jiyoye_essential_genes = jiyoye_essentiality_data.index
		raji_essential_genes = raji_essentiality_data.index

		gene_list = list(set(data.index))

		m = Mapper(gene_list,
				   input = 'symbol,uniprot,entrezgene,ensembl.gene,ensembl.protein',
				   output = 'symbol',
				   species = self.species)

		symbol_dict = m.hitDict['symbol']

		gene_list = list(data.index)
		raji_essentials = list()
		kbm7_essentials = list()
		k562_essentials = list()
		jiyoye_essentials = list()
		for gene in gene_list:
			try:
				sym = symbol_dict[gene]
				if sym.upper() in raji_essential_genes:
					raji_essentials.append(sym)
				else:
					raji_essentials.append('')
				if sym.upper() in kbm7_essential_genes:
					kbm7_essentials.append(sym)
				else:
					kbm7_essentials.append('')
				if sym.upper() in k562_essential_genes:
					k562_essentials.append(sym)
				else:
					k562_essentials.append('')
				if sym.upper() in jiyoye_essential_genes:
					jiyoye_essentials.append(sym)
				else:
					jiyoye_essentials.append('')
			except:
				raji_essentials.append('')
				kbm7_essentials.append('')
				k562_essentials.append('')
				jiyoye_essentials.append('')
		data['raji_essential'] = pd.Series(raji_essentials, index = data.index)
		data['kbm7_essential'] = pd.Series(kbm7_essentials, index = data.index)
		data['k562_essential'] = pd.Series(k562_essentials, index = data.index)
		data['jiyoye_essential'] = pd.Series(jiyoye_essentials, index = data.index)
		return data

	@staticmethod
	def load_string_data(string_folder, species):
		print('load string-data')
		file_name = species + "_STRING_geneName_per_protein_allInteractingProteins_dict.json"
		with open(string_folder + file_name) as json_data:
			string_dict_all = json.load(json_data)

		file_name = species + "_STRING_only500_geneName_per_protein_allInteractingProteins_dict.json"
		with open(string_folder + file_name) as json_data:
			string_dict_500 = json.load(json_data)

		file_name = species + "_STRING_only700_geneName_per_protein_allInteractingProteins_dict.json"
		with open(string_folder + file_name) as json_data:
			string_dict_700 = json.load(json_data)

		return string_dict_all, string_dict_500, string_dict_700

	@staticmethod
	def map_string_interactions(data, folder, species):
		string_dict_all, string_dict_500, string_dict_700 = step1_mapping.load_string_data(folder, species)

		gene_list = list(set(data.index))

		m = Mapper(gene_list,
				   input = 'symbol,uniprot,entrezgene,ensembl.gene,ensembl.protein',
				   output = 'symbol')

		symbol_dict = m.hitDict['symbol']

		gene_list = list(data.index)
		string_all_list = list()
		string_500_list = list()
		string_700_list = list()
		num_string_all_list = list()
		num_string_500_list = list()
		num_string_700_list = list()
		for gene in gene_list:
			try:
				sym = symbol_dict[gene]
				try:
					interactors = string_dict_all[sym]
					string_all_list.append(';'.join(interactors))
					num_string_all_list.append(len(interactors))
				except:
					string_all_list.append('')
					num_string_all_list.append('')
				try:
					interactors = string_dict_500[sym]
					string_500_list.append(';'.join(interactors))
					num_string_500_list.append(len(interactors))
				except:
					string_500_list.append('')
					num_string_500_list.append('')
				try:
					interactors = string_dict_700[sym]
					string_700_list.append(';'.join(interactors))
					num_string_700_list.append(len(interactors))
				except:
					string_700_list.append('')
					num_string_700_list.append('')
			except:
				string_all_list.append('')		
				string_500_list.append('')		
				string_700_list.append('')		
				num_string_all_list.append('')		
				num_string_500_list.append('')		
				num_string_700_list.append('')
		data['STRING_all'] = pd.Series(string_all_list, index = data.index)
		data['STRING_500'] = pd.Series(string_500_list, index = data.index)
		data['STRING_700'] = pd.Series(string_700_list, index = data.index)
		data['STRING_all_num'] = pd.Series(num_string_all_list, index = data.index)
		data['STRING_500_num'] = pd.Series(num_string_500_list, index = data.index)
		data['STRING_700_num'] = pd.Series(num_string_700_list, index = data.index)
		return data

	@staticmethod
	def map_chromosome_interactions(data, species):
		gene_list = list(set(data.index))

		if species == 'human':

			m = Mapper(gene_list,
					   input = 'symbol,uniprot,entrezgene,ensembl.gene,ensembl.protein',
				   	   output = 'genomic_pos_hg19.chr',
				   	   species = species)
			hit_dict = m.hitDict['genomic_pos_hg19']

		elif species == 'mouse':
			m = Mapper(gene_list,
					   input = 'symbol,uniprot,entrezgene,ensembl.gene,ensembl.protein',
				   	   output = 'genomic_pos_mm9.chr',
				   	   species = species)
			hit_dict = m.hitDict['genomic_pos_mm9']

		gene_list = list(data.index)
		chrom_list = list()
		for gene in gene_list:
			try:
				chrom_list.append(hit_dict[gene]['chr'])
			except:
				chrom_list.append('')
		data['chromosome'] = pd.Series(chrom_list, index = data.index)
		return data

	@staticmethod
	def export_dataframe(data, output_folder, output_fileName):
		data.to_csv(output_folder + output_fileName + '.tsv.gz', sep = '\t', compression = 'gzip')

class step1:
	@staticmethod
	def execute(**kwargs):
		folder = kwargs.get('folder','PATH')
		species  =  kwargs.get('species','human')#mouse available
		
		#THE OUTCOMMENTED STEPS HAVE BEEN EXECUTED ALREADY
		#step1_string_preparation.execute(folder = folder, species = species)
		#tcga_ovarian_cancer_formatting_prep(folder)
		#tcga_breast_cancer_formatting_prep(folder)

		fileName = "df_original_gygiData_tableS1_founderStrains.txt"
		step1.mapping_data(folder, fileName, 'Protein ID', 'mouse')
		fileName = "df_original_gygiData_tableS2_DOmice_RNA.txt"
		step1.mapping_data(folder, fileName, 'Gene ID','mouse')
		fileName = "df_original_gygiData_tableS3_DOmice_protein.txt"
		step1.mapping_data(folder, fileName, 'Protein ID', 'mouse')
		fileName =  "df_original_wuData_Analysis.txt"
		step1.mapping_data(folder, fileName, 'ENSPP', 'human')
		fileName = "df_original_battleProteinAnalysis.txt"
		step1.mapping_data(folder, fileName, 'ENSGG', 'human')
		fileName = "df_original_battleRiboAnalysis.txt"
		step1.mapping_data(folder, fileName, 'ENSGG', 'human')
		fileName = "df_original_battleRNAAnalysis.txt"
		step1.mapping_data(folder, fileName, 'ENSGG', 'human')
		fileName  =  "df_original_primateProteinData_Analysis.txt"
		step1.mapping_data(folder, fileName, 'ENSGG', 'human')
		fileName  =  "df_original_primateRNAData_Analysis.txt"
		step1.mapping_data(folder, fileName, 'ENSGG', 'human')
		fileName  =  "df_original_log2_tiannanData_Analysis.txt"		
		step1.mapping_data(folder, fileName, 'uniprot', 'human')
		fileName  =  "df_original_bxdMouse_dataset_proteomics.txt"
		step1.mapping_data(folder, fileName, 'GeneSymbol', 'mouse')
		fileName  =  "df_original_colorectal_cancer_dataset.txt"
		step1.mapping_data(folder, fileName, 'Gene name', 'human')
		fileName  =  "df_original_TCGA_Breast_ALL_Proteome.itraq.tsv"
		step1.mapping_data(folder, fileName, 'Gene name', 'human')
		fileName  =  "df_original_TCGA_Ovarian_ALL_Proteome.itraq.tsv"
		step1.mapping_data(folder, fileName, ' ', 'human')
		step1.mapping_geiger_HumanCelltypes_protein(folder)

	@staticmethod
	def mapping_data(folder, filename, proteinID, species):
		fileData  =  DataFrameAnalyzer.getFile(folder, fileName)
		outputFileName  =  "_".join(fileName.split("_")[1:]).split(".")[0]

		info1  =  step1_wider_mapping(fileData,
									  proteinID,
									  outputFileName = outputFileName,
									  species = species)

		info2  =  step1_module_mapping.execute(info1.mappedData,
											   folder = folder, 
											   species = species, 
											   output_folder = folder,
											   output_fileName = outputFileName)

	@staticmethod
	def mapping_geiger_HumanCelltypes_protein(folder):
		fileName = "df_original_log2_mann_11cell_Analysis.txt"
		fileData = DataFrameAnalyzer.getFile(folder, fileName)
		uniprotList = list(fileData["uniprot"])
		idxList = list()
		for uni in uniprotList:
			for item in uni.split(";"):
				if item.find("-") == -1:
					idxList.append(item)
					break
		fileData["uniprot"] = pd.Series(idxList,index = fileData.index)
		proteinID = "uniprot"
		outputFileName = "_".join(fileName.split("_")[1:]).split(".")[0]

		info1 = step1_wider_mapping(fileData,
									proteinID,
									outputFileName = outputFileName, 
									species = 'human')

		info2 = step1_module_mapping.execute(info1.mappedData,
											 folder = folder,
											 species = 'human', 
											 output_folder = folder,
											 output_fileName = outputFileName)

if __name__ == "__main__":
	## EXECUTE STEP 1
	step1.execute(folder = sys.argv[1])


