
# short-term projection model

run_projection <- function(report = report, par_est = par_est, )

#selectivity
n_sel_blocks = 5
sel_trawl_blocl = 3
sel_trap_block = 4

#growth block
GTM_block = c(1,2)

#natural mortality
M_year = 33
#M = report$Natural_Mortality[report$N_year*(report$Time_step-1)+M_year,]
M=0.4
SFrac = 1/3

n_rep = 10000
n_proj_year = 5

F_trawl_proj = rep(1.0,n_proj_year)
F_trap_proj  = rep(0.2,n_proj_year)

R_proj = c(500,500,500,500,500)

report <- read.admb(filename)
par_est <- std_read("nslsap01.std")

recr_lbin <- report$Recrut_Sizebin_ratio
recr_lbin_n <- length(recr_lbin[recr_lbin>0])
lbin_n <- length(report$Size_bins)-1
lbin <- report$Size_bins[1:lbin_n]+0.5

PropFemale = report$Maturity_length[report$N_year,]
weight = report$Weight_length[report$N_year,]

#GTM

GTM_1 = report$Growth[((GTM_block[1]-1)*lbin_n+1):((GTM_block[1]-1)*lbin_n+lbin_n),]
GTM_2 = report$Growth[((GTM_block[2]-1)*lbin_n+1):((GTM_block[2]-1)*lbin_n+lbin_n),]

sel_trawl_pars = c(par_est$est[(sel_trawl_blocl-1)*2+1][[1]],par_est$est[(sel_trawl_blocl-1)*2+2][[1]])
sel_trap_pars = c(par_est$est[(sel_trap_block-1)*2+1][[1]],par_est$est[(sel_trap_block-1)*2+2][[1]])
sel_trawl_pars_sd = c(par_est$std[(sel_trawl_blocl-1)*2+1][[1]],par_est$std[(sel_trawl_blocl-1)*2+2][[1]])
sel_trap_pars_sd = c(par_est$std[(sel_trap_block-1)*2+1][[1]],par_est$std[(sel_trap_block-1)*2+2][[1]])

N1_proj = matrix(NA,nrow=lbin_n,ncol=2)
N1_proj[(recr_lbin_n+1):lbin_n,1] = par_est$est$MvNASPJD_1
N1_proj[(recr_lbin_n+1):lbin_n,2] = par_est$std$MvNASPJD_1
N1_proj_temp = c()


SSB = matrix(NA,nrow=n_rep,ncol=n_proj_year)
C_trawl_tot=matrix(NA,nrow=n_rep,ncol=n_proj_year)
C_trap_tot=matrix(NA,nrow=n_rep,ncol=n_proj_year)
N_proj_3d=array(NA,dim=c(n_rep,lbin_n,n_proj_year))

for (y in 1:n_proj_year){
  for(i in 1:n_rep)
  {
    # define recruitment
    R_at_size_proj = R_proj[y]*recr_lbin
    
    # randomly draw selectivety parameter
    sel_trawl_pars_t=c(rnorm(1,sel_trawl_pars[1],sel_trawl_pars_sd[1]),rnorm(1,sel_trawl_pars[2],sel_trawl_pars_sd[2]))
    sel_trap_pars_t=c(rnorm(1,sel_trap_pars[1],sel_trap_pars_sd[1]),rnorm(1,sel_trap_pars[2],sel_trap_pars_sd[2]))
    sel_trawl = 1/(1+exp(-sel_trawl_pars_t[2]*(lbin-sel_trawl_pars_t[1])));sel_trawl=sel_trawl/max(sel_trawl)
    sel_trap  = 1/(1+exp(-sel_trap_pars_t[2]*(lbin-sel_trap_pars_t[1])));sel_trap=sel_trap/max(sel_trap)
    # define F  
    F_trawl = sel_trawl*F_trawl_proj[y]
    F_trap = sel_trap*F_trap_proj[y]
    
    # randomly draw N1_proj
    if (y==1){
      N1_proj_temp = R_at_size_proj
      for (ii in (recr_lbin_n+1):lbin_n)
        N1_proj_temp[ii] = rnorm(1,N1_proj[ii,1],N1_proj[ii,2])
      N_proj_3d[i,,y] = N1_proj_temp
    }else{
      N_proj_3d[i,,y] = (((N_proj_3d[sample(1:n_rep,1,FALSE),,y-1] * exp(-sel_trawl*F_trawl_proj[y-1]-sel_trap*F_trap_proj[y-1]-M))%*%GTM_1) *exp(-M))%*%GTM_2 + R_at_size_proj
      #N1_proj_temp = R_at_size_proj+N1_proj_temp
    }
    
    # population dynamic
    C_trawl         = N_proj_3d[i,,y]*(1-exp(-F_trawl-F_trap-M))*(F_trawl/(F_trawl+F_trap+M))
    C_trap          = N_proj_3d[i,,y]*(1-exp(-F_trawl-F_trap-M))*(F_trap/(F_trawl+F_trap+M))
    
    SSB[i,y]          = sum(N_proj_3d[i,,y]*exp(-SFrac*(-F_trawl-F_trap-M))*PropFemale*weight)
    C_trawl_tot[i,y]  = sum(C_trawl*weight)
    C_trap_tot[i,y]   = sum(C_trap*weight)
  
  }
}

#visulization 
library("fanplot")
par(mfrow=c(2,1))
par(mar = c(4, 4.5, 2, 0.5))
Year=seq(report$BeginYear,report$EndYear,1)
if(report$SexAtSizeLamda>0)
{SSB_estimates = report$Spawning_stock_Biomass
}else
  SSB_estimates=report$Spawning_stock_Biomass_input

plot(Year,SSB_estimates,xlab="Year",ylab="Spawning biomass (mt)",type='l',xlim=c(report$BeginYear-1,report$EndYear+n_proj_year))
fan(data = SSB, start = report$EndYear+1,rlab=F)

plot(NULL,xlab="Year",ylab="Spawning biomass (mt)",type='l',xlim=c(report$EndYear,report$EndYear+n_proj_year),ylim=c(500,2400))
fan(data = SSB, start = report$EndYear+1, ln=c(5, 50, 95))

Table1 = array(NA, dim=c(n_proj_year,3), dimnames=list(c(as.character(c(seq(1,n_proj_year,1)))),c("SSB","Catch_trawl","Catch_trap")) )

Table1[,1] = apply(SSB,2,median)
Table1[,2] = apply(C_trawl_tot,2,median)
Table1[,3] = apply(C_trap_tot,2,median)

print(Table1)













