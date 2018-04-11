
setwd("/Users/jiecao/Desktop/Umaine_model/Updated_runs/base_case")
source("reptoRlist.r");source("ADMB_read.r");source("run_projection.r")
filename="NSLSAP01"
report<-read.admb(filename);par_est <- std_read("nslsap01.std")
source("PlotFuncs.r")

pdf("plots.pdf",width=8,height=6)

PlotWL() # plot length-weight relationships
PlotML() # plot maturity
Plot_GM_tb(1,10) # visualize growth matrices (user defines the number of year-specific time-blocks for growth and how many year to grow)
Plot_sel_f() # plot fleet selectivity
Plot_f() # plot fishing mortality
PlotTC() # plot total catch
Plot_f_c() # plot fishery-specific catch
PlotCC(1,1) # plot fishery-specific catch composition; specify fishery and season
PlotCC(2,1)
PlotCC(3,1)
PlotCCA(3) # plot aggregated fishery catch composition
PlotSelS(2000,1) # plot survey selectivity; specify year and survey
PlotSelS(2000,2)
#PlotSelS(2000,3)
PlotSLC() # plot survey size composition
PlotSLCA() # plot aggregated survey size composition
PlotSI(3) # plot survey index;n= # of indices
PlotM(2000)
PlotR()
PlotSSB()
PlotAbun()
PlotSpB(1)
PlotSexComp(1,0.5)
PlotFfit()
PlotLfifty()
PlotER()

#PlotStockStatus()
#PlotRoSSB(1985,2002,2006,4)
#PlotRo(1984,2012,2016,2)
#PlotRo_scale(1984,2008,2013,2)
#Getlike()
dev.off()

# projection
run_projection(report = report, par_est = par_est, select_id = c(5,3,4), M = 0.25,
			    N_proj = 5, N_rep = 10000, F_trawl_proj = c(0,0,0,0,0), F_trap_proj = c(0,0,0,0,0), 					R_proj = c(1000,1000,1000,1000,1000), growth_id = c(1,2))














