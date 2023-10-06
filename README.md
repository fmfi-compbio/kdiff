# README
```sh
git clone --recursive git@github.com:ldenti/mfcnv.git
cd mfcnv
make
```

```
cd example
kmc -k21 -fq @case.lst case .
kmc -k21 -fq @control.lst control .

../main reference.fa control case 100 > bins.depth 2> pos.depth
```
