#!/bin/sh

module load nco

var_list="Times,XLONG,XLAT,HFX,LH,FIRA,FSA,GPP,SWDOWN,FVEG,RSSUN,RSSHA,W,T,P,T2V,T2B,TAH,TV,TH2,T2,PH,PHB,HGT,U,V,PSFC,U10,V10"
#var_list="Times,XLONG,XLAT,W,T,P,T2V,T2B,TAH,TV,TH2,T2"

yatirdir="/global/cscratch1/sd/twhilton/WRFv4.1_Experiments/WRFv4.1_yatir_2015Aug_yatirparams_NCEPFNL_NOAHMP_ndom3_wetsoil/WRFV4/run/yatir_2015Aug_yatirparams_NCEPFNL_NOAHMP_ndom3_wetsoil"
#yatirdir="/global/cscratch1/sd/twhilton/WRFv4.1_Experiments/WRFv4.1_yatir_2015Aug_yatirparams_Z050_NCEPFNL_NOAHMP_ndom3/WRFV4/run/yatir_2015Aug_yatirparams_NCEPFNL_NOAHMP_ndom3"
ctldir="/global/cscratch1/sd/twhilton/WRFv4.1_Experiments/WRFv4.1_yatir_2015Aug_CTL_NCEPFNL_NOAHMP_ndom3/WRFV4/run/yatir_2015Aug_CTL_NCEPFNL_NOAHMP_ndom3"
outdir="/global/cscratch1/sd/twhilton/yatir_output_collected/wetsoil"

echo "processing Yatir run from $yatirdir to $outdir"
fname_yatir="yatir_run_d03_diag_TP_VWCx2.nc"
rm -fv $outdir/$fname_yatir
cd $yatirdir
echo "starting yatir `date`"
ncrcat -v $var_list metem_yatir_2015Aug_yatirparams_NCEPFNL_NOAHMP_ndom3_wetsoil_d03_2015-08-15_0[012]* $outdir/$fname_yatir
echo "done yatir `date`"

echo "processing control run from $ctldir to $outdir"
fname_ctl="ctl_run_d03_diag_TP_VWCx2.nc"
rm -fv $outdir/$fname_ctl
cd $ctldir
echo "starting control `date`"
ncrcat -v $var_list metem_yatir_2015Aug_CTL_NCEPFNL_NOAHMP_ndom3_d03_2015-08-15_0[012]* $outdir/$fname_ctl
echo "done control `date`"
