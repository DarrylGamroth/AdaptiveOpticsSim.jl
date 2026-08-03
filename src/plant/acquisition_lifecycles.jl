"""
Common cold declaration boundary for scheduled acquisition lifecycles.

Detector lifecycles and intentional direct-measurement lifecycles share this
boundary without pretending that a direct WFS measurement is a detector
readout.
"""
abstract type AbstractAcquisitionLifecycleDefinition end

"""
    AbstractPreparedAcquisitionLifecycle

Qualified-public prepared-owner interface for one scheduled acquisition.
Implementations provide `begin_exposure!`, `complete_readout!`, exact-target
validation through `_require_exact_acquisition_lifecycle_target`, and
`structural_resource_fact` reporting. A prepared lifecycle binds one exact
definition and product/resource owner; its matching mutable state is
single-writer and non-reentrant. Every transition rejects a foreign binding or
invalid transition with a structured Plant or detector-acquisition error before
mutation. Independent prepared/state pairs may execute independently.
"""
abstract type AbstractPreparedAcquisitionLifecycle end

"""Common mutable single-writer state boundary for an acquisition lifecycle."""
abstract type AbstractAcquisitionLifecycleState end
