# Generate graphs that compare results for different coefficients. 

rm(list=ls())
library(tictoc)
library(xtable)
library(doParallel)
n_cores = detectCores()
cl <- makeCluster(n_cores-1, setup_strategy = "sequential")
registerDoParallel(cl)


options("scipen"=100, "digits"=4)
tic()

d = 20
n = 1e05
m="m02"

#n.RData = 2
n.RData = 100
itenumber = 1000
n.cut = seq(25000, 100000, length.out=4)

dd = c(1:5)

coverage.rs.1 = matrix(0, nrow=n, ncol=1)
coverage.bm.1 = matrix(0, nrow=n, ncol=1)
coverage.pi.1 = matrix(0, nrow=n, ncol=1)

CI.length.rs.1 = matrix(0, nrow=n, ncol=1)
CI.length.bm.1 = matrix(0, nrow=n, ncol=1)
CI.length.pi.1 = matrix(0, nrow=n, ncol=1)

time.rs = matrix(0, nrow=n, ncol=1)
time.bm = matrix(0, nrow=n, ncol=1)
time.pi = matrix(0, nrow=n, ncol=1)


for (i in c(1:n.RData)){
  tic()
  cat("Proceeding......",i,'\n')
  load(paste0("logistic_d",d,"_",m,"_",i,".RData"))
  coverage.rs.1 = coverage.rs.1 + r.rs$covrate_ave
  coverage.bm.1 = coverage.bm.1 + r.bm$covrate_ave
  coverage.pi.1 = coverage.pi.1 + r.pi$covrate_ave
  
  CI.length.rs.1 = CI.length.rs.1 + r.rs$cilengthbar
  CI.length.bm.1 = CI.length.bm.1 + r.bm$cilengthbar
  CI.length.pi.1 = CI.length.pi.1 + r.pi$cilengthbar
  
  time.rs = time.rs + r.rs$runtime_cum
  time.bm = time.bm + r.bm$runtime_cum
  time.pi = time.pi + r.pi$runtime_cum
  
  gc()
  toc()
}

ave.coverage.rs = coverage.rs.1 / n.RData
ave.coverage.bm = coverage.bm.1 / n.RData
ave.coverage.pi = coverage.pi.1 / n.RData

coverage.rs.sd = sqrt( (ave.coverage.rs*(1-ave.coverage.rs)) / itenumber )
coverage.bm.sd = sqrt( (ave.coverage.bm*(1-ave.coverage.bm)) / itenumber )
coverage.pi.sd = sqrt( (ave.coverage.pi*(1-ave.coverage.pi)) / itenumber )

ave.CI.length.rs = CI.length.rs.1 / n.RData
ave.CI.length.bm = CI.length.bm.1 / n.RData
ave.CI.length.pi = CI.length.pi.1 / n.RData



ave.runtime.rs = time.rs / n.RData
ave.runtime.bm = time.bm / n.RData
ave.runtime.pi = time.pi / n.RData



# Parameters for graphs
col = c("darkred", "darkblue", "black")  # rs: red, bm: blue, pi: purple
lwd = 2
lty = c(1,2,4)

lr=1
alpha=0.505
c = 1
burn = 5
ak = as.integer(c*((1:10))**(2/(1-alpha)))
n.start = ak[burn] 



# Draw the first coefficient only: coverage rate
pdf(file = paste0("fig-coverage-d",d,"-",m,".pdf"), width = 8, height = 6) # The height of the plot in inches

plot(c(1:n), ave.coverage.rs, type="l", lty=1, col=col[1], lwd=lwd, ylim=c(0.6,1),
     main="Coverage Rates", xlab="Obs. Number", ylab="Coverage Rate",
     cex.lab=1.5, cex.axis=1.5, cex.main=2)

lines(ave.coverage.bm, type="l", lty=2, col=col[2], lwd=lwd )
lines(ave.coverage.pi, type="l", lty=4, col=col[3], lwd=lwd )
abline(h=0.95, lty=1, col="darkgreen", lwd=lwd )

legend(70000, 0.7, legend=c("Random Scale", "Batch-mean", "Plug-in"), col=col, lty=lty, cex=1.2, lwd=lwd)
dev.off()


# Draw the first coefficient only: CI length
pdf(file = paste0("fig-length-d",d,"-",m,".pdf"), width = 8, height = 6) # The height of the plot in inches

plot(c(1:n), ave.CI.length.rs, type="l", lty=1, col=col[1], lwd=lwd, ylim=c(0,0.1),
     main="CI Length", xlab="Obs. Number", ylab="CI Length",
     cex.lab=1.5, cex.axis=1.5, cex.main=2)

lines(ave.CI.length.bm, type="l", lty=2, col=col[2], lwd=lwd )
lines(ave.CI.length.pi, type="l", lty=4, col=col[3], lwd=lwd )

legend(70000, 0.1, legend=c("Random Scale", "Batch-mean", "Plug-in"), col=col, lty=lty, cex=1.2, lwd=lwd)
dev.off()

# Draw the first coefficient only: computation time
pdf(file = paste0("fig-time-d",d,"-",m,".pdf"), width = 8, height = 6) # The height of the plot in inches

plot(c(1:n), ave.runtime.rs, type="l", lty=1, col=col[1], lwd=lwd, ylim=c(0,60),
     main="Computation Time", xlab="Obs. Number", ylab="Time (sec.)",
     cex.lab=1.5, cex.axis=1.5, cex.main=2)

lines(ave.runtime.bm, type="l", lty=2, col=col[2], lwd=lwd )
lines(ave.runtime.pi, type="l", lty=4, col=col[3], lwd=lwd )

legend(0, 60, legend=c("Random Scale", "Batch-mean", "Plug-in"), col=col, lty=lty, cex=1.2, lwd=lwd)
dev.off()




table.name = paste0("table_d_",d)  

# Rows for coverage
coverage.rs = rbind(paste0(sprintf("%.3f",ave.coverage.rs[n.cut[1:4]])," (",sprintf("%.4f",coverage.rs.sd[n.cut[1:4]]),")") )
coverage.bm = rbind(paste0(sprintf("%.3f",ave.coverage.bm[n.cut[1:4]])," (",sprintf("%.4f",coverage.bm.sd[n.cut[1:4]]),")") )
coverage.pi = rbind(paste0(sprintf("%.3f",ave.coverage.pi[n.cut[1:4]])," (",sprintf("%.4f",coverage.pi.sd[n.cut[1:4]]),")") )
coverage = rbind(coverage.rs, coverage.bm, coverage.pi)

# Rows for length
#CI.length.rs = rbind(paste0(sprintf("%.3f",ave.CI.length.rs[n.cut[1:4]])," (",sprintf("%.4f",CI.length.rs.sd[n.cut[1:4]]),")") )
#CI.length.bm = rbind(paste0(sprintf("%.3f",ave.CI.length.bm[n.cut[1:4]])," (",sprintf("%.4f",CI.length.bm.sd[n.cut[1:4]]),")") )
#CI.length.pi = rbind(paste0(sprintf("%.3f",ave.CI.length.pi[n.cut[1:4]])," (",sprintf("%.4f",CI.length.pi.sd[n.cut[1:4]]),")") )
CI.length.rs = rbind(paste0(sprintf("%.3f",ave.CI.length.rs[n.cut[1:4]])) )
CI.length.bm = rbind(paste0(sprintf("%.3f",ave.CI.length.bm[n.cut[1:4]])) )
CI.length.pi = rbind(paste0(sprintf("%.3f",ave.CI.length.pi[n.cut[1:4]])) )
CI.length = rbind(CI.length.rs, CI.length.bm, CI.length.pi)

time = rbind(sprintf("%.1f",ave.runtime.rs[n.cut]), sprintf("%.1f",ave.runtime.bm[n.cut]), sprintf("%.1f",ave.runtime.pi[n.cut]))
table = rbind(coverage, CI.length, time)
print(xtable(table, caption = paste0("Design: ","lr=",lr,", alpha=",alpha)), type="latex", file=paste0(table.name,".tex"), include.rownames=FALSE)

toc()
stopCluster(cl)
