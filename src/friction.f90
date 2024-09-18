module friction
    !! Provides a collection of routines for modeling frictional behaviors
    !! of contacting bodies.
    use friction_core
    use friction_coulomb
    use friction_lugre
    use friction_maxwell
    use friction_gmsm
    use friction_stribeck
    use friction_modified_stribeck
    implicit none
    private
    public :: friction_model
    public :: friction_evaluation
    public :: friction_logical_query
    public :: friction_state_model
    public :: friction_model_to_array
    public :: friction_model_from_array
    public :: friction_integer_query
    public :: regression_statistics
    public :: coulomb_model
    public :: lugre_model
    public :: maxwell_model
    public :: generalized_maxwell_slip_model
    public :: stribeck_model
    public :: modified_stribeck_model

end module