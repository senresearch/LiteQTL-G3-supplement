using CUDA
# using CuArrays
using Base
using BenchmarkTools
using LinearAlgebra
using Dates
using Random

include("benchmark.jl")

const alf = 1.0
const bet = 0.0
const gb =  1073741824
const ns2sec = 1.0e-9

# gemm!(tA, tB, alpha, A, B, beta, C)

function cpu_run(a::Array, b::Array, c::Array)
    return LinearAlgebra.BLAS.gemm!('N','N',alf, a,b, bet, c)
end

# CuArrays.CUBLAS.gemm!('N','N',alpha,d_A,d_B,beta,d_C1) 

function gpu_run(a::Array, b::Array, c::Array)
    A = CuArray(a)
    B = CuArray(b)
    C = CuArray(c)
    CUDA.CUBLAS.gemm!('N', 'N', alf, A,B, bet, C)
    return collect(C)
end


# dt_now = Dates.format(Dates.now(), "yyyy-mm-ddTHH-MM-SS")
# host = gethostname()

timing_record_file = joinpath(@__DIR__, "..", "data", "gemm-timing", "dgemm-result.csv")

file = open(timing_record_file, "w")
thread_counts = 16 #[1,16,32]
for t in thread_counts
    LinearAlgebra.BLAS.set_num_threads(t)
    core_nums = ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
    println("\nnumber of threads used: $core_nums")
    write(file, "# number of cores/threads used in openblas: $core_nums\n");
    close(file)
    msizes = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576];
    nsizes = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576]; # range from 1k to 1m in log scale
    psizes = [16, 32, 64, 128, 512, 1024, 2048, 4096, 8192, 16384, 32768]; # ranges from 10 to 4k in log scale.

    for m in msizes
        for p in psizes
            for n in nsizes
                file = open(timing_record_file, "a")
                req_mem_gb = (m*n + n*p + m*p) * sizeof(Float64) / gb
                avail_mem_gb = CUDA.available_memory() / gb
                if (req_mem_gb > avail_mem_gb)
                    
                    println("Matrices are too big to fit in GPU memory. Skipping this configuration. M is $m, N is $n, P is $p");
                    write(file, "$m,$n,$p,$req_mem_gb,$avail_mem_gb,NA,NA,NA\n");
                    close(file);
                else
                    println("m = $m, n = $n, p: $p\n")
                    Random.seed!(123);

                    #generating double precision matrix
                    A = randn(m,n);
                    B = randn(n,p);
                    C = zeros(m,p);

                    println("CPU runtime:")
                    cpuelapsed = @elapsed cpu_run(A,B,C);
                    println("CPU took $(cpuelapsed) sec to run")

                    println("GPU runtime:")
                    gpuelapsed = @elapsed CUDA.@sync gpu_run(A,B,C);
                    println("GPU took $(gpuelapsed) s to run")


                    cpu_result = benchmark(10, cpu_run, A,B,C)
                    gpu_result = benchmark(10, gpu_run,A,B,C)
                    speedup = cpu_result[3]/gpu_result[3]
                    println(cpu_result)
                    println(gpu_result)

                    write(file, "$m,$n,$p,$req_mem_gb,$avail_mem_gb,$(cpu_result[3]),$(gpu_result[3]),$speedup\n");
                    close(file)
                end
            end
        end
    end
end
# file = open(timing_record_file, "a")
dt_finish = Dates.format(Dates.now(), "yyyy-mm-ddTHH-MM-SS")
write(file, "Finshied at $dt_finish");
close(file)
