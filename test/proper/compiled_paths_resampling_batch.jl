# Independent centered-DFT oracle for staged image resampling, including
# isolation of batch planes and a final partial group on changed-input replay.
function check_compiled_image_resampling_batch(array, target; captured=false)
    backend=AdaptiveOpticsSim.Backends
    context=backend._prepare_device_execution_context(target)
    backend._with_prepared_device_execution_context(context) do
        @testset "Image resampling batch $array" begin
            for T in (Float32,Float64), n in (5,6)
                side=n+2
                kernel=Complex{T}[cis(T(0.19i*j+0.03i^2-0.07j)) for i in 1:n,j in 1:n]
                plan=(;target,image_size=n,crop_start=1,n_anchor=1,
                    resample_kernel=array(kernel),
                    paths=map(i->(;path_index=i,anchor_index=1),(3,1,2)))
                closed=CP._compiled_paths_seal([array(zeros(Complex{T},n,n)) for _ in 1:3])
                only=CP._compiled_paths_seal([array(zeros(Complex{T},n,n)) for _ in 1:3])
                field=array(zeros(Complex{T},side,side))
                coefficients=backend._prepare_mul_coefficients(field)
                batch=CP._compiled_paths_prepare_image_resample(plan,closed,only,coefficients,Val(4))
                owner=(;plan,workspace=(;field))
                inputs=ntuple(_->array(zeros(Complex{T},side,side)),6)
                operation()=begin
                    for i in 1:6
                        backend.copyto_backend_async!(field,inputs[i])
                        slot=plan.paths[cld(i,2)].path_index
                        out=isodd(i) ? closed[slot] : only[slot]
                        CP._compiled_paths_stage_image!(batch,owner,i,out,true)
                    end
                end
                recorded=captured ? backend._capture_prepared_device_graph(operation,context) : nothing
                for state in 1:3
                    host=ntuple(6) do i
                        active=state==1 ? i==1 : state==2 ? i==6 : true
                        Complex{T}[active ? Complex{T}(sin(T(.17y+.23x+i)),cos(T(.13y-.11x+state))) : zero(Complex{T})
                            for y in 1:side,x in 1:side]
                    end
                    for i in 1:6
                        copyto!(inputs[i],array(host[i]))
                    end
                    fill!(batch.field,Complex{T}(NaN))
                    fill!(batch.temporary,Complex{T}(NaN))
                    recorded===nothing ? operation() : backend._launch_prepared_device_graph!(recorded,context)
                    backend._synchronize_prepared_device_execution_context!(context)
                    for i in 1:6
                        crop=host[i][2:n+1,2:n+1]
                        expected=(kernel*fftshift(fft(ifftshift(crop)))*transpose(kernel))/n^2
                        slot=plan.paths[cld(i,2)].path_index
                        actual=Array(isodd(i) ? closed[slot] : only[slot])
                        @test norm(actual-expected)/max(norm(expected),one(T))<=512eps(T)
                        (state==1 && i!=1 || state==2 && i!=6) && @test all(iszero,actual)
                    end
                    @test all(isfinite,Array(batch.field)) && all(isfinite,Array(batch.temporary))
                end
            end
        end
    end
end
