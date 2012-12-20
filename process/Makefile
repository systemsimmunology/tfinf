TFINF = ~/tfinf
PROC_DIR = $(TFINF)/process
OUT_DIR = $(TFINF)/data/process_output

vpath %.R $(PROC_DIR)
## Trying to get vpath to work has been a headache, particuraly for both PROC_DIR and OUT_DIR
## Try to fix later. See e.g. http://mad-scientist.net/make/vpath.html
## http://stackoverflow.com/questions/5178125/gnu-make-how-to-get-object-files-in-separate-subdirectory

all: $(addprefix $(OUT_DIR)/, tteMaps.RData TFcategories.RData expressed.scanned.ensembl.RData initializeBC.RData featureMatrix.mf.RData ESmaxcubeMgrid.RData SingleFeatures.RData PairedFeatures.RData sigPairedSitesBC.RData pdna.curated.RData pdnaModels.RData scaled.mus.objects.RData boost.vec.RData collection.1reg.RData collection.2regs.RData models.rmsf.RData modvec.e.RData)

$(OUT_DIR)/tteMaps.RData: expressedTFandMats.R
	R --no-restore --no-save < expressedTFandMats.R >& expressedTFandMats.log

$(OUT_DIR)/TFcategories.RData: TFcategories.R $(OUT_DIR)/tteMaps.RData
	R --no-restore --no-save < TFcategories.R >& TFcategories.log

$(OUT_DIR)/expressed.scanned.ensembl.RData $(OUT_DIR)/Parsed.Scan.eset.RData: expressed.scanned.ensembl.R
	R --no-restore-data --slave  < expressed.scanned.ensembl.R

$(OUT_DIR)/initializeBC.RData: initializeBC.R $(OUT_DIR)/Parsed.Scan.eset.RData
	R --no-restore --no-save < initializeBC.R

$(OUT_DIR)/featureMatrix.mf.RData : $(OUT_DIR)/initializeBC.RData find.pairs_allMouseBC.R
	R --no-save < find.pairs_allMouseBC.R > find.pairs_allMouseBC.logfile

$(OUT_DIR)/ESmaxcubeMgrid.RData: $(OUT_DIR)/initializeBC.RData randEScube.R
	R --no-save < randEScube.R > randEScube.logfile

$(OUT_DIR)/SingleFeatures.RData:  GetFeatureListSingles.R $(OUT_DIR)/featureMatrix.mf.RData $(OUT_DIR)/ESmaxcubeMgrid.RData
	R --no-save < GetFeatureListSingles.R > GetFeatureListSingles.logfile

$(OUT_DIR)/PairedFeatures.RData: GetFeatureListPairs.R $(OUT_DIR)/featureMatrix.RData $(OUT_DIR)/ESmaxcubeMgrid.RData
	R --no-save < GetFeatureListPairs.R >  GetFeatureListPairs.logfile

$(OUT_DIR)/sigPairedSitesBC.RData: metacollectionWrapperBC.R $(OUT_DIR)/SingleFeatures.RData $(OUT_DIR)/PairedFeatures.RData $(OUT_DIR)/featureMatrix.mf.RData $(OUT_DIR)/featureMatrix.RData
	R --no-save < metacollectionWrapperBC.R > metacollectionWrapperBC.logfile

$(OUT_DIR)/pdna.curated.RData: curatedPDNAints.R
	R --no-restore --no-save < curatedPDNAints.R > curatedPDNAints.logfile

$(OUT_DIR)/pdnaModels.RData: createTargsAndsCandsBC.R $(OUT_DIR)/sigPairedSitesBC.RData $(OUT_DIR)/pdna.curated.RData
	R --no-restore --no-save < createTargsAndsCandsBC.R > createTargsAndsCandsBC.logfile

$(OUT_DIR)/scaled.mus.objects.RData: insertNucLoc.R
	R --no-restore --no-save < insertNucLoc.R

$(OUT_DIR)/boost.vec.RData: createBoostVec.R ## Requires boost.vec.refit.RData in AUX_DIR
	R --no-restore --no-save < createBoostVec.R

$(OUT_DIR)/collection.1reg.RData: fixedTargsRandomTFsOneReg.R $(OUT_DIR)/boost.vec.RData
	R --no-save < fixedTargsRandomTFsOneReg.R >& fixedTargsRandomTFsOneReg.log

$(OUT_DIR)/collection.2regs.RData: fixedTargsRandomTFsTwoRegs.R $(OUT_DIR)/boost.vec.RData
	R --no-save < fixedTargsRandomTFsTwoRegs.R >& fixedTargsRandomTFsTwoRegs.log

$(OUT_DIR)/models.rmsf.RData $(OUT_DIR)/models.ode.RData: calcModelsWrapperFinal.R
	R --no-restore --no-save < calcModelsWrapperFinal.R >& calcModelsWrapperFinal.logfile

$(OUT_DIR)/modvec.e.RData: combineModelClasses.R $(OUT_DIR)/models.rmsf.RData
	R --no-restore --no-save < combineModelClasses.R
