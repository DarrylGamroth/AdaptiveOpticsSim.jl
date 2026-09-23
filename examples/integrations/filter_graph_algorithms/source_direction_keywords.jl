# The frozen integration environment pins an AOS revision with the historical
# source-direction keyword. Keep that spelling at the fixture boundary only.
@inline source_direction_keywords(::Val{:frozen}, separation_arcsec, position_angle_deg) =
    (coordinates=(separation_arcsec, position_angle_deg),)

@inline source_direction_keywords(::Val{:current}, separation_arcsec, position_angle_deg) =
    (; separation_arcsec, position_angle_deg)
