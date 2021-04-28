using CSV 
using Statistics


liteqtl_timing = CSV.read(timing_file, comment="#",missingstring="NA")

lqtl_cpu = liteqtl_timing[liteqtl_timing[:device] .== "CPU",:]
# first ten data entry is for full matrix timing 
lqtl_fullmat_cpu = lqtl_cpu[1:10, :]
# the rest data entry is for filter result timing
lqtl_filtered_cpu = lqtl_cpu[11:end, :]

lqtl_gpu = liteqtl_timing[liteqtl_timing[:device] .== "GPU",:]
# first ten data entry is for full matrix timing 
lqtl_fullmat_gpu = lqtl_gpu[1:10, :]
# the rest data entry is for filter result timing
lqtl_filtered_gpu = lqtl_gpu[11:end, :]

println("elapsed runtime (- pval_calc_time) =  (data transfer) + (compute time) + (result reorg time)")
println("Full Matrix: CPU $(median(lqtl_fullmat_cpu.elapsed_total)) = $(median(lqtl_fullmat_cpu.data_transfer_time)) + $(median(lqtl_fullmat_cpu.compute_time)) + $(median(lqtl_fullmat_cpu.result_reorg_time))")
println("Full Matrix: GPU $(median(lqtl_fullmat_gpu.elapsed_total)) = $(median(lqtl_fullmat_gpu.data_transfer_time)) + $(median(lqtl_fullmat_gpu.compute_time)) + $(median(lqtl_fullmat_gpu.result_reorg_time))") 
println("Max Only: CPU $(median(lqtl_filtered_cpu.elapsed_total)) = $(median(lqtl_filtered_cpu.data_transfer_time)) + $(median(lqtl_filtered_cpu.compute_time)) + $(median(lqtl_filtered_cpu.result_reorg_time))")
println("Max Only: GPU $(median(lqtl_filtered_gpu.elapsed_total)) = $(median(lqtl_filtered_gpu.data_transfer_time)) + $(median(lqtl_filtered_gpu.compute_time)) + $(median(lqtl_filtered_gpu.result_reorg_time))")

# tensorqtl_timing = CSV.read("tensorqtl_timing_report.txt",comment="#")
# #get all cpu timing.
# tqtl_cpu = tensorqtl_timing[tensorqtl_timing[:device] .== "cpu",:]
# # first ten data entry is for full matrix timing 
# tqtl_fullmat_cpu = tqtl_cpu[1:10, :]
# # the rest data entry is for filter result timing
# tqtl_filtered_cpu = tqtl_cpu[11:end, :]

# tqtl_gpu = tensorqtl_timing[tensorqtl_timing[:device] .== "cuda",:]
# # first ten data entry is for full matrix timing 
# tqtl_fullmat_gpu = tqtl_gpu[1:10, :]
# # the rest data entry is for filter result timing
# tqtl_filtered_gpu = tqtl_gpu[11:end, :]

# println("TensorQTL: ")
# println("elapsed runtime (- pval_calc_time) =  (data transfer) + (compute time) + (result reorg time)")
# println("Full Matrix: CPU: $(median(tqtl_fullmat_cpu.elapsed_total) - median(tqtl_fullmat_cpu.pval_time)), $(median(tqtl_fullmat_cpu.data_transfer_time)), $(median(tqtl_fullmat_cpu.compute_time)),$(median(tqtl_fullmat_cpu.result_reorg_time))")
# println("Full Matrix: GPU: $(median(tqtl_fullmat_gpu.elapsed_total) - median(tqtl_fullmat_gpu.pval_time)), $(median(tqtl_fullmat_gpu.data_transfer_time)),$(median(tqtl_fullmat_gpu.compute_time)),$(median(tqtl_fullmat_gpu.result_reorg_time))") 
# println("Filtered: CPU: $(median(tqtl_filtered_cpu.elapsed_total) - median(tqtl_filtered_cpu.pval_time)), $(median(tqtl_filtered_cpu.data_transfer_time)),$(median(tqtl_filtered_cpu.compute_time)),$(median(tqtl_filtered_cpu.result_reorg_time))")
# println("Filtered: GPU: $(median(tqtl_filtered_gpu.elapsed_total) - median(tqtl_filtered_gpu.pval_time)), $(median(tqtl_filtered_gpu.data_transfer_time)),$(median(tqtl_filtered_gpu.compute_time)),$(median(tqtl_filtered_gpu.result_reorg_time))")
