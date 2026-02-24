"""Named physical dimension constants for the floorplan.

All values in feet unless noted. Inches converted via / 12.0.
"""

# Shell construction (constant regardless of wall thickness)
SHELL_THICKNESS = 2.0 / 12.0         # 2" concrete shell

# Wall thicknesses (feet)
WALL_OUTER = 8.0 / 12.0           # 8" outer wall (adjustable: 8"-12")
assert 8.0 / 12.0 <= WALL_OUTER <= 12.0 / 12.0, "WALL_OUTER must be 8\"-12\""
WALL_EXTRA = WALL_OUTER - 8.0 / 12.0     # wall extra beyond 8" baseline
WALL_6IN = 6.0 / 12.0             # 6" interior wall (IW1, IW2, IW2s)
WALL_4IN = 4.0 / 12.0             # 4" interior wall (IW3, IW4)
WALL_3IN = 3.0 / 12.0             # 3" interior wall (IW5)

# Appliance dimensions (feet)
APPLIANCE_WIDTH = 35.0 / 12.0     # 35" washer/dryer width
APPLIANCE_DEPTH = 30.0 / 12.0     # 30" washer/dryer depth
APPLIANCE_OFFSET_FROM_W2 = 6.0 / 12.0   # 6" from W2, CW-normal to W2-W5
APPLIANCE_OFFSET_FROM_W1 = 4.0 / 12.0   # 4" from W1, along W2-W5 direction
APPLIANCE_GAP = 1.0 / 12.0        # 1" gap between dryer and washer

# Counter
COUNTER_DEPTH = 24.0 / 12.0       # 2' E-W
COUNTER_LENGTH = 6.0              # 6' N-S
COUNTER_GAP = 36.0 / 12.0         # 3' east of dryer

# Rooms
IW1_OFFSET_FROM_W1 = 12.0 + 8.0/12.0    # 12'8" from W1, along W2-W5 direction
IW1_OFFSET_FROM_W2 = 6.5           # 6'6" from W2, CW-normal to W2-W5
# IW2 segments (breakup of former single IW2 into lower, oblique, shower)
IW2_DIST_W2W5 = 6.5                  # 6'6" from W2-W5 inner wall (true perpendicular)
IW2_LENGTH = 42.0 / 12.0             # 42" north from IW1 north face
IW2S_W2REF_OFFSET = 6.5              # offset from virtual W2 ref → ~5'4" true from W2-W5
IW2S_LENGTH = 72.0 / 12.0            # 72" (6') south from W6-W7 inner wall
IW2O_THICKNESS = 6.0 / 12.0          # 6" thick perpendicular to midline
IW4_GAP_IW11 = 30.0 / 12.0           # 30" from IW11 east face to IW4 west face
RO1_OFFSET_FROM_IW2 = 116.0 / 12.0   # 9'8" from IW2 east face, CW-normal to W2-W5

# Bed
BED_WIDTH = 76.0 / 12.0           # 76" king bed
BED_LENGTH = 94.0 / 12.0          # 94" (incl. frame)
BED_WALL_GAP = 2.0 / 12.0         # 2" bed-to-outer-wall inward gap

# Dresser
DRESSER_WIDTH = 34.0 / 12.0       # 34" E-W
DRESSER_DEPTH = 19.0 / 12.0       # 19" N-S
DRESSER_GAP_IW15 = 2.0 / 12.0     # 2" west of IW15
DRESSER_GAP_IW1 = 1.0 / 12.0      # 1" south of IW1

# Water heater
WH_RADIUS = 14.0 / 12.0           # 14" radius (28" diameter)

# Toilet (plan view)
TOILET_WIDTH = 15.0 / 12.0        # 15"
TOILET_TANK_DEPTH = 8.0 / 12.0    # 8" tank

# Sink (plan view, ellipse semi-axes)
SINK_RX = 12.0 / 12.0             # 12" E-W half-width (24" total)
SINK_RY = 9.0 / 12.0              # 9" N-S half-depth (18" total)

# Kitchen appliances (feet)
KITCHEN_SINK_WIDTH = 45.0 / 12.0   # 45" E-W
KITCHEN_SINK_DEPTH = 24.0 / 12.0   # 24" N-S
DW_WIDTH = 28.0 / 12.0             # 28" dishwasher E-W
DW_DEPTH = 27.0 / 12.0             # 27" dishwasher N-S
STOVE_WIDTH = 30.0 / 12.0          # 30" stove E-W
STOVE_DEPTH = 27.0 / 12.0          # 27" stove N-S
FRIDGE_SIZE = 36.0 / 12.0          # 36" fridge (square)
MINIK_FRIDGE_W = 23.375 / 12.0    # 23-3/8" IKEA BERGSNAS width (E-W)
MINIK_FRIDGE_D = 24.75 / 12.0     # 24-3/4" IKEA BERGSNAS depth (N-S)
KITCHEN_GAP = 0.75 / 12.0          # 3/4" gap between kitchen appliances
KITCHEN_CTR_LENGTH = 72.0 / 12.0   # 72" kitchen counter E-W along IW1 north
KITCHEN_CTR_DEPTH = 30.0 / 12.0    # 30" kitchen counter depth N-S
NORTH_CTR_LENGTH = 36.0 / 12.0     # 36" north wall counter E-W
NORTH_CTR_DEPTH = 30.0 / 12.0      # 30" north wall counter depth N-S
EAST_CTR_LENGTH = 30.0 / 12.0      # 30" east counter E-W along W9-W10
EAST_CTR_DEPTH = 42.0 / 12.0       # 42" east counter depth N-S

# IW1 rough opening
RO1_OFFSET_FROM_IW9 = 76.0 / 12.0    # 76" from IW9 east face
IW1_RO_WIDTH = 38.0 / 12.0        # 38" opening width E-W

# IW2o rough opening (RO4 centered on IW2o oblique segment)
IW2_RO_WIDTH = 38.0 / 12.0        # 38" opening width along IW2o

# IW3 (perpendicular to W20-W0, 4" thick)
IW3_LENGTH = 80.0 / 12.0          # 80" (6'8") length into structure
IW3_OFFSET_IW9 = 30.0 / 12.0     # 30" from IW9 W face along W20-W0
IW3_DIST_W2W5 = 102.0 / 12.0     # 8'6" from W2-W5 to IW3 west face

# IW7 (parallel to W20-W0, 4" thick, between IW3 and IW9)

# IW9 (perpendicular to W20-W0, 4" thick)
IW9_LENGTH = 80.0 / 12.0            # 80" (6'8") IW9 length, same as IW3
IW9_OFFSET_O10 = 6.0 / 12.0         # 6" past O10 along W20-W0

# IW11 (4" thick, N-S)
IW9_IW11_GAP = 12.0                  # 12' from IW9 east face to IW11 west face

# IW12 (4" thick, perpendicular to IW11)
IW12_OFFSET_IW11 = 6.0              # 6' from IW11 SW to IW12 base
IW12_SHORTEN = 4.0 / 12.0           # 4" IW12 west-end setback
IW9_RO_WIDTH = 62.0 / 12.0          # 62" opening width along IW9 (RO7)

# IW16 rough opening
IW16_RO_WIDTH = 38.0 / 12.0         # 38" opening width N-S

# IW4 rough opening
IW4_RO_WIDTH = 38.0 / 12.0        # 38" opening width N-S

# IW11 rough opening (RO6)
IW11_RO_WIDTH = 62.0 / 12.0       # 62" opening width along IW11

# IW6 partition
IW6_THICKNESS = 1.0 / 12.0        # 1" partition
IW6_OFFSET_FROM_W6 = 5.5                 # 5'6" from W6, CW-normal to W6-W7

# IW6 rough opening
IW6_RO_OFFSET_W = 3.0 / 12.0      # 3" west of IW2 west face
IW6_RO_WIDTH = 38.0 / 12.0        # 38" opening width E-W

# Outer-wall openings (numbered CW around outline)
# O1 (F2-F5 east wall, centered at RO3 normal projection)
O1_WIDTH = 19.0 / 12.0             # 19" opening height
# O2 (F2-F5, centered at RO4 northing center)
O2_WIDTH = 19.0 / 12.0             # 19" opening width
# O3 (F2-F5, 4" from F5)
O3_GAP_F5 = 4.0 / 12.0             # 4" from F5 along F5-F2 line
O3_WIDTH = 32.0 / 12.0             # 32" opening width
O3_DOOR_WIDTH = 30.0 / 12.0        # 30" door in O3
# O4 (F6-F7, relative to IW2 west face)
O4_HALF_WIDTH = 4.5 / 12.0         # 4.5" half-width (9" total)
O4_OFFSET_FROM_IW2 = 11.0 / 12.0     # 11" from IW2
# O5 (F9-F10)
O5_OFFSET_FROM_IW2 = 120.0 / 12.0       # 10' from IW2, CW-normal to W2-W5
O5_WIDTH = 68.0 / 12.0            # 5'8" opening width
# O6 (F9-F10)
O6_WIDTH = 44.0 / 12.0             # 44" opening width
O6_GAP_F10 = 6.0 / 12.0            # 6" from O6 east edge to F10
# F10 easting: 15'2" east of nominal F9
F10_OFFSET_FROM_F9 = 182.0 / 12.0     # 15'2" from F9, along W9-W10
O6_DOOR_WIDTH = 42.0 / 12.0        # 42" door, centered in opening
RO1_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO1
RO2_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO2
RO3_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO3
RO4_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO4
RO5_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO5
RO6_DOOR_WIDTH = 30.0 / 12.0       # 30" door leaf in RO6 (double door, 2×30")
RO7_DOOR_WIDTH = 30.0 / 12.0       # 30" door leaf in RO7 (double door, 2×30")
# Shell construction: 2" shell / gap / 2" shell
AIR_GAP = WALL_OUTER - 2 * SHELL_THICKNESS  # air gap between shells
OPENING_INSIDE_RADIUS = 10.0 / 304.8  # 10mm inside corner radius at openings
# Door jamb block thickness = wall - 2*(opening_inside_radius + shell_thickness)
DOOR_FLAT_FACE = WALL_OUTER - 2 * (OPENING_INSIDE_RADIUS + SHELL_THICKNESS)
# F8-F9 inner wall turn radius (W-face = inner face of inner shell)
F8F9_INNER_TURN_R = OPENING_INSIDE_RADIUS + SHELL_THICKNESS  # ~2.56" (10mm + 2")
# O7 (F12-F13 diagonal wall)
O7_NW_GAP = 2.0                    # 2' from F12 to NW end
O7_HALF_WIDTH = 36.0 / 12.0        # 36" half-width (72" total opening)
# O8 (F14-F15)
O8_HALF_WIDTH = 11.5 / 12.0        # 11.5" half-width
# O9, O10, O11 (F20-F1 south wall chain)
O9_HALF_WIDTH = 11.5 / 12.0        # 11.5" half-width
O10_HALF_WIDTH = 11.5 / 12.0       # 11.5" half-width
O11_HALF_WIDTH = 9.5 / 12.0        # 9.5" half-width
O9_OFFSET_IW11 = 6.0 / 12.0       # 6" IW11 SW to O9 SE along F20-F1
O9_O10_WALL = 86.0 / 12.0         # 86" solid wall between O9 NW and O10 SE
O10_O11_WALL = 72.0 / 12.0        # 72" solid wall between O10 NW and O11 SE
BED_GAP_O9 = 4.0 / 12.0           # 4" from O9 NW to bed SE along W20-W1

# IW5 partition
IW5_OFFSET_FROM_IW1 = 30.0 / 12.0        # 30" from IW1, CW-normal to W9-W10

# Outline geometry constraints
CORNER_NE_R = 10.0 / 12.0 + WALL_EXTRA   # R_a1: NE corner (10" at 8" wall)
F11AB_TARGET = 1.0                 # 1'0" target F11a-F11b distance
IW1_OFFSET_FROM_W9 = 11.0              # 11'0" from W9, CW-normal to W9-W10
IW8_OFFSET_FROM_W20W1 = 12.0          # 12'0" perpendicular to W20-W1

# Jamb and gap constants
JAMB_WIDTH = 1.0 / 12.0           # 1" jamb width (rough openings)
STD_GAP = 2.0 / 12.0              # 2" standard gap (furniture/appliance spacing)
KITCHEN_APPL_GAP = 3.0 / 12.0     # 3" gap (kitchen appliance spacing/setback)

# Furniture dimensions
LOVESEAT_WIDTH = 35.0 / 12.0      # 35" loveseat short side
LOVESEAT_LENGTH = 65.0 / 12.0     # 65" loveseat long side
LOVESEAT_OFFSET_IW4 = 4.522741716102669  # feet along IW4 outward from IW4/IW1 corner
LOVESEAT_OFFSET_IW1 = 5.609454568566491  # feet along IW1 outward from IW4/IW1 corner
DESK_WIDTH = 60.0 / 12.0          # 60" desk (along wall)
DESK_DEPTH = 30.0 / 12.0          # 30" desk (perpendicular to wall)
DESK_CHAIR_WIDTH = 27.0 / 12.0   # 27" desk chair E-W
DESK_CHAIR_DEPTH = 24.0 / 12.0   # 24" desk chair N-S
DESK_CHAIR_GAP = 12.0 / 12.0     # 12" gap between desk and chair
CHAIR_WIDTH = 32.0 / 12.0         # 32" chair E-W
CHAIR_DEPTH = 37.0 / 12.0         # 37" chair N-S
CHAIR_CORNER_R = 3.0 / 12.0       # 3" rounded corner radius
OTTOMAN_SIZE = 29.0 / 12.0        # 29" square ottoman
ET_RADIUS_CM = 25.0               # 25 cm endtable radius
SOFA_WIDTH = 97.2 / 12.0          # 97.2" sofa E-W
SOFA_DEPTH = 24.6 / 12.0          # 24.6" sofa N-S
ICE_WIDTH = 17.7 / 12.0           # 17.7" ice maker E-W
ICE_DEPTH = 15.8 / 12.0           # 15.8" ice maker N-S
ROCKER_WIDTH = 26.75 / 12.0       # 26.75" POANG rocking chair E-W
ROCKER_DEPTH = 37.0 / 12.0        # 37" POANG rocking chair N-S
ROCKER_CORNER_R = 3.0 / 12.0      # 3" rounded corner radius

# Roof overhang
ROOF_OVERHANG = 6.0 / 12.0        # 6" roof overhang beyond wall face
