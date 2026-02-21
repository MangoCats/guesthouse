"""Named physical dimension constants for the floorplan.

All values in feet unless noted. Inches converted via / 12.0.
"""

# Shell construction (constant regardless of wall thickness)
SHELL_THICKNESS = 2.0 / 12.0         # 2" concrete shell

# Wall thicknesses (feet)
WALL_OUTER = 8.0 / 12.0           # 8" outer wall (adjustable: 8"-12")
assert 8.0 / 12.0 <= WALL_OUTER <= 12.0 / 12.0, "WALL_OUTER must be 8\"-12\""
_WE = WALL_OUTER - 8.0 / 12.0     # wall extra beyond 8" baseline
WALL_6IN = 6.0 / 12.0             # 6" interior wall (IW1, IW2)
WALL_4IN = 4.0 / 12.0             # 4" interior wall (IW3, IW4)
WALL_3IN = 3.0 / 12.0             # 3" interior wall (IW7, IW8)

# Appliance dimensions (feet)
APPLIANCE_WIDTH = 35.0 / 12.0     # 35" washer/dryer width
APPLIANCE_DEPTH = 30.0 / 12.0     # 30" washer/dryer depth
APPLIANCE_OFFSET_E = 6.0 / 12.0   # 6" from west wall to dryer
APPLIANCE_OFFSET_N = 4.0 / 12.0   # 4" from south wall to dryer
APPLIANCE_GAP = 1.0 / 12.0        # 1" gap between dryer and washer

# Counter
COUNTER_DEPTH = 24.0 / 12.0       # 2' E-W
COUNTER_LENGTH = 70.0 / 12.0      # 5'10" N-S
COUNTER_NW_RADIUS = 9.0 / 12.0    # 9" rounded corner
COUNTER_GAP = 36.0 / 12.0         # 3' east of dryer

# Rooms
BEDROOM_WIDTH = 138.0 / 12.0      # 11'6" E-W
CLOSET_WIDTH = 30.0 / 12.0        # 30" closet depth (closet 1)
CLOSET2_WIDTH = 28.0 / 12.0       # 28" closet 2 depth (east)
CLOSET1_HEIGHT = 6.0              # 6' closet 1 N-S
IW1_OFFSET_N = 12.0 + 8.0/12.0    # 12'8" IW1 south face above W1
IW1_WEST_OFFSET_E = 6.5           # 6'6" IW1 west end east of W2
IW2_OFFSET_E = 6.5                # 6'6" IW2 west face east of W2
IW4_OFFSET_E_IW2 = 224.0 / 12.0   # 18'8" IW4 west face east of IW2 east face
RO1_OFFSET_E_IW2 = 116.0 / 12.0   # 9'8" RO1 west edge east of IW2 east face
WALL_SOUTH_N = 2.0 / 12.0         # 2" south end of bedroom walls

# Bed
BED_WIDTH = 76.0 / 12.0           # 76" king bed
BED_LENGTH = 94.0 / 12.0          # 94" (incl. frame)
BED_OFFSET_N = 2.0 / 12.0         # 2" from south wall

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
EAST_CTR_RADIUS = 12.0 / 12.0      # 12" south corner radius

# IW1 rough opening
RO1_OFFSET_W_IW4 = 64.0 / 12.0    # 64" west of IW4 west face (legacy)
RO1_OFFSET_E_IW9 = 76.0 / 12.0    # 76" east of IW9 east face
IW1_RO_WIDTH = 38.0 / 12.0        # 38" opening width E-W

# IW2 rough opening
IW2_RO_OFFSET_S = 9.0 / 12.0      # 9" south of IW6 south face
IW2_RO_WIDTH = 38.0 / 12.0        # 38" opening width N-S

# IW3 (perpendicular to W20-W0, 4" thick)
IW3_LENGTH = 80.0 / 12.0          # 80" (6'8") length into structure
IW3_OFFSET_IW9 = 30.0 / 12.0     # 30" from IW9 W face along W20-W0

# IW7 (parallel to W20-W0, 4" thick, between IW3 and IW9)

# IW9 (perpendicular to W20-W0, 4" thick)
IW9_LENGTH = 80.0 / 12.0            # 80" (6'8") IW9 length, same as IW3
IW9_OFFSET_O10 = 6.0 / 12.0         # 6" past O10 along W20-W0

# IW16 rough opening
IW16_RO_WIDTH = 38.0 / 12.0         # 38" opening width N-S

# IW4 rough opening
IW4_RO_WIDTH = 38.0 / 12.0        # 38" opening width N-S

# IW6 partition
IW6_THICKNESS = 1.0 / 12.0        # 1" partition
IW6_OFFSET_N = 5.5                 # 5'6" south of F6-F7 south face

# IW6 rough opening
IW6_RO_OFFSET_W = 3.0 / 12.0      # 3" west of IW2 west face
IW6_RO_WIDTH = 38.0 / 12.0        # 38" opening width E-W

# Outer-wall openings (numbered CW around outline)
# O1 (F2-F3, lower)
O1_OFFSET_S = 99.0 / 12.0          # 99" south of F3 to north edge
O1_WIDTH = 19.0 / 12.0             # 19" opening height
# O2 (F4-F5, centered at RO4 northing center)
O2_WIDTH = 19.0 / 12.0             # 19" opening width
# O3 (F4-F5, 4" from F5)
O3_GAP_F5 = 4.0 / 12.0             # 4" from F5 along F5-F4 line
O3_WIDTH = 32.0 / 12.0             # 32" opening width
O3_DOOR_WIDTH = 30.0 / 12.0        # 30" door in O3
# O4 (F6-F7, relative to IW2 west face)
O4_HALF_WIDTH = 4.5 / 12.0         # 4.5" half-width (9" total)
O4_OFFSET_W_IW2 = 11.0 / 12.0     # 11" west of IW2 west face
# O5 (F9-F10)
O5_E_FROM_IW2 = 120.0 / 12.0       # 10' from IW2 east face to O5 east edge
O5_WIDTH = 68.0 / 12.0            # 5'8" opening width
# O6 (F9-F10)
O6_WIDTH = 44.0 / 12.0             # 44" opening width
O6_GAP_F10 = 6.0 / 12.0            # 6" from O6 east edge to F10
# F10 easting: 15'2" east of nominal F9
F10_OFFSET_E_F9 = 182.0 / 12.0     # 15'2" from nominal F9 to F10
O6_DOOR_WIDTH = 42.0 / 12.0        # 42" door, centered in opening
RO1_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO1
RO2_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO2
RO3_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO3
RO4_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO4
RO5_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO5
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
IW5_OFFSET_N = 30.0 / 12.0        # 30" south of IW1 south face

# Outline geometry constraints
CORNER_NE_R = 10.0 / 12.0 + _WE   # R_a1: corner arc (10" at 8" wall)
CORNER_NW_R = 28.0 / 12.0 + _WE   # R_a5: NW corner (28" at 8" wall)
UPPER_E_R = 28.0 / 12.0 + _WE     # R_a7: upper east (28" at 8" wall)
SMALL_ARC_R = 2.0 / 12.0          # R_a8: 2" transition
ARC_180_R = 28.0 / 12.0 + _WE     # R_a11: 180-degree arc (28" at 8" wall)
FLAT_SEG_11 = 16.0 / 12.0         # 16" straight segment F11a-F11b
ARC_F3_R = 6.064223608163559 + _WE  # R_a3: ~72.8" at 8" wall
F6_EAST_ADJ = 6.0 / 12.0           # 6" F6 east adjustment
F6_HEIGHT = 27.0 + 2*_WE  # F6-F7 north of F1 (27'0" at 8"; +2*_WE because F1 moves south and F6 north)
NW_SHIFT = 1.0                    # C5 1' east shift
IW1_DIST_FROM_NORTH = 11.0              # 11'0" IW1 north face to north inner wall
IW8_OFFSET_N_IW1 = 19.0 / 12.0    # 19" IW8 north face above IW1 north face
F3_GAP_N_IW8 = 4.0 / 12.0         # 4" F3 north of IW8 north face (legacy)
ARC_F3_SWEEP = 10.0                # 10° F3-F4 arc sweep (F5-F6 = 90° - this)
F14_OFFSET_N_IW1 = 2.0 / 12.0     # 2" F14 north of IW1 north face
F14_F15_SEG = 8.5                  # 8'6" segment
ARC_F13_R = 28.0 / 12.0 + _WE        # R_a13: 28" at 8" wall
ARC_F13_R_BASELINE = 5.627004870830987 + _WE  # C11a anchor baseline (original R_a13)
F13_EXIT_BRG = 345.0              # 345-degree exit bearing
SOUTH_WALL_N = -6.0 / 12.0 - _WE  # south face wall northing (-6" at 8" wall)
PIX_PI5_TARGET_BRG = 60.0         # 60-degree target bearing
F15_OFFSET_E = 9.0 + 3.0/12.0 + _WE  # F15 east of iw8_e (9'3" at 8" wall)
ARC_F17_SWEEP = 30.0               # 30° sweep for F17-F18 arc (CW)
F16_F17_MIN = 5.0                  # minimum 5' F16-F17 segment length
F18_OFFSET_E = 2.0 / 12.0         # 2" min F18 east of IW4 east face
F19_OFFSET_E = -10.0 / 12.0       # F19 10" west of IW4 east face
ARC_F19_R = 18.888718471469218 + _WE  # R_a19: ~226.7" at 8" wall

# Jamb and gap constants
JAMB_WIDTH = 1.0 / 12.0           # 1" jamb width (rough openings)
STD_GAP = 2.0 / 12.0              # 2" standard gap (furniture/appliance spacing)
KITCHEN_APPL_GAP = 3.0 / 12.0     # 3" gap (kitchen appliance spacing/setback)

# Furniture dimensions
LOVESEAT_WIDTH = 35.0 / 12.0      # 35" loveseat short side
LOVESEAT_LENGTH = 65.0 / 12.0     # 65" loveseat long side
LOVESEAT_ANGLE_DEG = 15.0         # 15° CCW rotation
LOVESEAT_NW_E = 22.310591617230667  # NW corner easting (fixed position)
LOVESEAT_NW_N = 18.94278790189982   # NW corner northing (fixed position)
DESK_WIDTH = 60.0 / 12.0          # 60" desk (along wall)
DESK_DEPTH = 30.0 / 12.0          # 30" desk (perpendicular to wall)
DESK_CHAIR_WIDTH = 27.0 / 12.0   # 27" desk chair E-W
DESK_CHAIR_DEPTH = 24.0 / 12.0   # 24" desk chair N-S
DESK_CHAIR_GAP = 12.0 / 12.0     # 12" gap between desk and chair
CHAIR_WIDTH = 32.0 / 12.0         # 32" chair E-W
CHAIR_DEPTH = 37.0 / 12.0         # 37" chair N-S
CHAIR_CORNER_R = 3.0 / 12.0       # 3" rounded corner radius
CHAIR_ANGLE_DEG = 30.0            # 30° CW rotation
OTTOMAN_SIZE = 29.0 / 12.0        # 29" square ottoman
ET_RADIUS_CM = 25.0               # 25 cm endtable radius
SOFA_WIDTH = 97.2 / 12.0          # 97.2" sofa E-W
SOFA_DEPTH = 24.6 / 12.0          # 24.6" sofa N-S
SHELVES_WIDTH = 36.0 / 12.0       # 36" shelves E-W
SHELVES_DEPTH = 15.0 / 12.0       # 15" shelves N-S
ICE_WIDTH = 17.7 / 12.0           # 17.7" ice maker E-W
ICE_DEPTH = 15.8 / 12.0           # 15.8" ice maker N-S
ROCKER_WIDTH = 26.75 / 12.0       # 26.75" POANG rocking chair E-W
ROCKER_DEPTH = 37.0 / 12.0        # 37" POANG rocking chair N-S
ROCKER_CORNER_R = 3.0 / 12.0      # 3" rounded corner radius

# Roof overhang
ROOF_OVERHANG = 6.0 / 12.0        # 6" roof overhang beyond wall face
