"""
Common cold declaration boundary for scheduled acquisition lifecycles.

Detector lifecycles and intentional direct-measurement lifecycles share this
boundary without pretending that a direct WFS measurement is a detector
readout.
"""
abstract type AbstractAcquisitionLifecycleDefinition end

"""Common prepared boundary for one scheduled acquisition lifecycle."""
abstract type AbstractPreparedAcquisitionLifecycle end

"""Common mutable single-writer state boundary for an acquisition lifecycle."""
abstract type AbstractAcquisitionLifecycleState end
