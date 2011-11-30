#!/bin/bash
R --no-save < find.pairs_allMouseBC.R > find.pairs_allMouseBC.logfile
R --no-save < randEScube.R > randEScube.logfile
R --no-save < GetFeatureListSingles.R > GetFeatureListSingles.logfile
R --no-save < GetFeatureListPairs.R >  GetFeatureListPairs.logfile
R --no-save < metacollectionWrapperBC.R > metacollectionWrapperBC.logfile
cd ../R
R --no-restore --no-save < createTargsAndsCandsBC.R > createTargsAndsCandsBC.logfile
