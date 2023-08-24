# Code for Tau3Mu-related analyses

## Getting started

```shell
cmsrel CMSSW_12_4_11
cd CMSSW_12_4_11/src
cmsenv
git cms-init
```

## checkout 
```
git clone git@github.com:BasChiara/Tau3MuNANO.git ./PhysicsTools/Tau3MuNANO
cd PhysicsTools/Tau3MuNANO
```

## make sure we use a consistent tag
```
git fetch origin
git checkout -b myBranch origin/myBranch
```

## add your own fork as a remote. Skip if you dont have one
```
git remote add crovelli git@github.com:crovelli/Tau3MuNANO.git
git fetch crovelli
```

## compile
cd ../../
scramv1 b