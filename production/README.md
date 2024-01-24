# instructions to submit on CRAB
source CRAB3 and activate proxy
`
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init --voms cms --valid 168:00
`
submit on dataset listed in ```data_samples_2023BC.yml```
`
python3 submit_on_crab.py -y data_samples_2023BC.yml
`
