# Code for Tau3Mu-related analyses

## Getting started

```shell
cmsrel CMSSW_13_0_13
cd CMSSW_13_0_13/src
cmsenv
git cms-init
```

## checkout 
```
git clone git@github.com:BasChiara/Tau3muNANO.git ./PhysicsTools/Tau3muNANO
cd PhysicsTools/Tau3muNANO
```

## make sure we use a consistent tag
```
git fetch origin
git checkout -b myBranch origin/myBranch
```

## add your own fork as a remote. Skip if you dont have one
```
git remote add crovelli git@github.com:crovelli/Tau3muNANO.git
git fetch crovelli
#git checkout -b devChiara crovelli/master
```

## compile
```
cd ../../
scramv1 b
```
