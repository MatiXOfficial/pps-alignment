# pps-alignment

## Installation
```
scram project CMSSW_11_2_0_pre2
cd CMSSW_11_2_0_pre2/src
cmsenv
git clone https://github.com/MatiXOfficial/pps-alignment.git .
scram b -j 8
```
## Run test
```
cd CalibPPS/Alignment/test
cmsRun test_PPSAlignment_cfg.py
```
## Event Setup example
```
cd CalibPPS/ESProducers/test
cmsRun test_PPSAlignmentConfigESSource_cfg.py
```