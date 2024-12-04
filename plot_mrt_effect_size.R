plot_mrt_effect_size <- function(mrt_effect_size_answer) {
  
  plot(x=mrt_effect_size_answer$time, 
       y=mrt_effect_size_answer$estimate, 
       lwd=2, 
       type="l", 
       lty="solid",
       xlab="Decision Time",
       ylab=expression("d"),
       main="Effect Size Estimate")
  lines(x=mrt_effect_size_answer$time, 
        y=mrt_effect_size_answer$lower,
        lwd=2,
        lty="dashed")
  lines(x=mrt_effect_size_answer$time, 
        y=mrt_effect_size_answer$upper,
        lwd=2,
        lty="dashed")
  
} 