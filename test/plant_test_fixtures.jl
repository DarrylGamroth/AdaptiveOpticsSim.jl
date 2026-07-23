function captured_plant_preparation_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function assert_plant_preparation_error(f, component::Symbol,
    reason::Symbol)
    error = captured_plant_preparation_error(f)
    @test error isa PlantPreparationError
    if error isa PlantPreparationError
        @test error.component === component
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

function prepared_selection_execution_allocations(selection,
    epoch::AtmosphereEpoch)
    execute_acquisition_selection!(selection, epoch)
    return @allocated execute_acquisition_selection!(selection, epoch)
end
