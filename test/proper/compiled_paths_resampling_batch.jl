# Independent centered-DFT oracle for staged image resampling, including
# physical finalization, isolation of batch planes and changed-input replay.
function check_compiled_image_resampling_batch(array, target; captured=false)
    backend=AdaptiveOpticsSim.Backends
    context=backend._prepare_device_execution_context(target)
    backend._with_prepared_device_execution_context(context) do
        @testset "Image resampling batch $array" begin
            for T in (Float32,Float64), n in (5,6), finalization_kind in (:none,:normalization,:post)
                side=n+2
                kernel=Complex{T}[cis(T(0.19i*j+0.03i^2-0.07j)) for i in 1:n,j in 1:n]
                post=finalization_kind===:post ?
                    Complex{T}[Complex{T}(0.7+0.01y+0.02x,0.1sin(0.2x-0.3y)) for y in 1:side,x in 1:side] :
                    zeros(Complex{T},0,0)
                normalization=side^2
                finalization=finalization_kind===:none ? nothing : (normalization,array(post))
                plan=(;target,image_size=n,crop_start=1,n_anchor=1,
                    resample_kernel=array(kernel),
                    paths=map(i->(;path_index=i,anchor_index=1),(3,1,2)))
                closed=CP._compiled_paths_seal([array(zeros(Complex{T},n,n)) for _ in 1:3])
                only=CP._compiled_paths_seal([array(zeros(Complex{T},n,n)) for _ in 1:3])
                field=array(zeros(Complex{T},side,side))
                finalized=similar(field)
                coefficients=backend._prepare_mul_coefficients(field)
                batch=CP._compiled_paths_prepare_image_resample(plan,closed,only,coefficients,Val(4))
                owner=(;plan,workspace=(;field))
                inputs=ntuple(_->array(zeros(Complex{T},side,side)),6)
                operation()=begin
                    for i in 1:6
                        backend.copyto_backend_async!(field,inputs[i])
                        slot=plan.paths[cld(i,2)].path_index
                        out=isodd(i) ? closed[slot] : only[slot]
                        CP._compiled_paths_stage_image!(batch,owner,i,out,true,finalization)
                    end
                end
                operation()
                backend._synchronize_prepared_device_execution_context!(context)
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
                    # Compare the pack exactly against the original full-field
                    # backend arithmetic: divide before multiplying the post.
                    if finalization===nothing
                        copyto!(finalized,inputs[1])
                    elseif isempty(post)
                        finalized .= inputs[1] ./ normalization
                    else
                        resident_post=last(finalization)
                        @. finalized=(inputs[1]/normalization)*resident_post
                    end
                    fill!(batch.field,Complex{T}(NaN))
                    CP._compiled_paths_pack_resample_kernel!(
                        CP.KernelAbstractions.get_backend(batch.field),(16,16))(
                        batch.field,inputs[1],plan.crop_start,2,2,finalization;ndrange=(n,n))
                    backend._synchronize_prepared_device_execution_context!(context)
                    packed=Array(batch.field)
                    expected_pack=ifftshift(Array(finalized)[2:n+1,2:n+1])
                    if T===Float64 && finalization_kind===:post && target!=backend.HostComputeDevice()
                        # Verified GPU contraction differs only in the complex
                        # post's imaginary sum: fma(b,c,a*d) vs fma(a,d,b*c).
                        # Scale by the operands, including at cancellations;
                        # retain a zero bound wherever either operand is zero.
                        normalized_crop=ifftshift(host[1][2:n+1,2:n+1]./normalization)
                        post_crop=ifftshift(post[2:n+1,2:n+1])
                        bound=4eps(T).*abs.(normalized_crop).*abs.(post_crop)
                        @test all(abs.(packed[:,:,2].-expected_pack).<=bound)
                    else
                        @test packed[:,:,2] == expected_pack
                    end
                    @test all(isnan,packed[:,:,1])
                    @test all(iszero,packed[:,:,3:4])
                    fill!(batch.field,Complex{T}(NaN))
                    fill!(batch.temporary,Complex{T}(NaN))
                    recorded===nothing ? operation() : backend._launch_prepared_device_graph!(recorded,context)
                    backend._synchronize_prepared_device_execution_context!(context)
                    for i in 1:6
                        physical=copy(host[i])
                        if finalization!==nothing
                            physical ./= normalization
                            isempty(post) || (physical .*= post)
                        end
                        crop=physical[2:n+1,2:n+1]
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
