CA.@calculon_algorithm ModalOPDExpansionF32 begin
    label = "modal-opd-expansion-f32"

    @configuration ModalOPDExpansionF32Configuration begin
        pupil_rows::Int => graph_rebuild
        pupil_columns::Int => graph_rebuild
        mode_count::Int => graph_rebuild
        coefficients_schema::String => graph_rebuild
        opd_schema::String => graph_rebuild
        basis_schema::String => graph_rebuild
        pupil_support_schema::String => graph_rebuild
    end

    properties = :none

    ports(configuration, plan, workspace) = (
        input(
            :coefficients,
            Float32,
            (configuration.mode_count,),
            configuration.coefficients_schema,
        ),
        output(
            :opd,
            Float32,
            (configuration.pupil_rows, configuration.pupil_columns),
            configuration.opd_schema;
            metadata_from=:coefficients,
        ),
        parameter(
            :basis,
            Float32,
            (
                configuration.pupil_rows,
                configuration.pupil_columns,
                configuration.mode_count,
            ),
            configuration.basis_schema;
            replace=(plan, values) ->
                Calibration._replace_modal_opd_basis(plan, values),
        ),
        parameter(
            :pupil_support,
            Bool,
            (configuration.pupil_rows, configuration.pupil_columns),
            configuration.pupil_support_schema;
            replace=(plan, values) ->
                Calibration._replace_modal_opd_pupil_support(plan, values),
        ),
    )

    prepare(configuration) = (
        Calibration.ModalOPDExpansionPlan(
            zeros(
                Float32,
                configuration.pupil_rows,
                configuration.pupil_columns,
                configuration.mode_count,
            ),
            fill(false, configuration.pupil_rows, configuration.pupil_columns),
        ),
        CA.NoWorkspace(),
    )

    process(plan, workspace) =
        Calibration.combine_basis!(opd, plan, coefficients)
end

public ModalOPDExpansionF32
public ModalOPDExpansionF32Configuration
