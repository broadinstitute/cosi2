#!/idi/sabeti-data/ilya/sabetilab/ilya/usr/bin/bash -e

N=60
#M=/idi/sabeti-data/ilya/nsvn/Operations/Ilya_Operations/sim/cosi21_gc_msHOTmodel/cosi/smp.cosiParams
#M=/idi/sabeti-data/ilya/nsvn/Operations/Ilya_Operations/sim/cosi21_gc_msHOTmodel/cosi/0_simple.cosiParams
M=/idi/sabeti-data/ilya/nsvn/Operations/Ilya_Operations/sim/cosi21_gc_msHOTmodel/cosi/cmplx.cosiParams

myvar=0
while [ $myvar -ne $N ]
do
	 ./coalescent -p $M
	 mv trajnew.tsv trajnew$myvar.tsv
    myvar=$(( $myvar + 1 ))
done

pushd /idi/sabeti-data/ilya/nsvn/Operations/Ilya_Operations/sim/cosi/cosi_befopt/cosi
myvar=0
while [ $myvar -ne $N ]
do
	 ./coalescent -p $M
	 mv trajold.tsv trajold$myvar.tsv
    myvar=$(( $myvar + 1 ))
done

cd /idi/sabeti-data/ilya/nsvn/Temp
python cmptraj.py $N
popd


