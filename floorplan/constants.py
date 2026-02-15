"""Named physical dimension constants for the floorplan.

All values in feet unless noted. Inches converted via / 12.0.
"""

# Wall thicknesses (feet)
WALL_OUTER = 8.0 / 12.0           # 8" outer wall
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
COUNTER_LENGTH = 72.0 / 12.0      # 6' N-S
COUNTER_NW_RADIUS = 9.0 / 12.0    # 9" rounded corner
COUNTER_GAP = 36.0 / 12.0         # 3' east of dryer

# Rooms
BEDROOM_WIDTH = 140.0 / 12.0      # 11'8" E-W
CLOSET_WIDTH = 30.0 / 12.0        # 30" closet depth
CLOSET1_HEIGHT = 6.0              # 6' closet 1 N-S
IW1_OFFSET_N = 11.5               # 11'6" IW1 south face above W0
IW2_OFFSET_E = 6.5                # 6'6" IW2 west face east of W1
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
KITCHEN_GAP = 0.75 / 12.0          # 3/4" gap between kitchen appliances
KITCHEN_CTR_LENGTH = 72.0 / 12.0   # 72" kitchen counter E-W along IW1 north
KITCHEN_CTR_DEPTH = 30.0 / 12.0    # 30" kitchen counter depth N-S
NORTH_CTR_LENGTH = 36.0 / 12.0     # 36" north wall counter E-W
NORTH_CTR_DEPTH = 30.0 / 12.0      # 30" north wall counter depth N-S
EAST_CTR_LENGTH = 30.0 / 12.0      # 30" east counter E-W along W9-W10
EAST_CTR_DEPTH = 42.0 / 12.0       # 42" east counter depth N-S
EAST_CTR_RADIUS = 12.0 / 12.0      # 12" south corner radius

# IW1 rough opening
IW1_RO_OFFSET_E = 9.0 / 12.0      # 9" east of fridge east side
IW1_RO_WIDTH = 38.0 / 12.0        # 38" opening width E-W

# IW2 rough opening
IW2_RO_OFFSET_S = 6.0 / 12.0      # 6" south of IW6 south face
IW2_RO_WIDTH = 38.0 / 12.0        # 38" opening width N-S

# IW3 rough opening
IW3_RO_OFFSET_N = 2.0 / 12.0      # 2" north of IW7 north face
IW3_RO_WIDTH = 38.0 / 12.0        # 38" opening width N-S

# IW4 rough opening
IW4_RO_WIDTH = 38.0 / 12.0        # 38" opening width N-S

# IW6 partition
IW6_THICKNESS = 1.0 / 12.0        # 1" partition
IW6_OFFSET_N = 5.5                 # 5'6" south of F6-F7 south face

# IW6 rough opening
IW6_RO_OFFSET_W = 3.0 / 12.0      # 3" west of IW2 west face
IW6_RO_WIDTH = 38.0 / 12.0        # 38" opening width E-W

# Outer-wall openings (numbered CW around outline)
# O1 (F1-F2, lower)
O1_OFFSET_S = 99.0 / 12.0          # 99" south of F2 to north edge
O1_WIDTH = 25.0 / 12.0             # 25" opening height
# O2 (F1-F2, upper)
O2_OFFSET_S = 4.0 / 12.0           # 4" south of F2 to north edge
O2_WIDTH = 25.0 / 12.0             # 25" opening height
# O3 (F4-F5, centered)
O3_HALF_WIDTH = 16.0 / 12.0        # 16" half-width
# O4 (F6-F7, centered chimney opening)
O4_HALF_WIDTH = 4.5 / 12.0         # 4.5" half-width (9" total)
# O5 (F9-F10)
O5_E_FROM_F7 = 108.0 / 12.0        # 9' from F7 easting to O5 east edge
O5_WIDTH = 68.0 / 12.0            # 5'8" opening width
# O6 (F9-F10)
O6_E_FROM_F9 = 194.0 / 12.0        # 16'2" from F9 to O6 east edge
O6_WIDTH = 44.0 / 12.0             # 44" opening width
F10_O6_CLEARANCE = 4.0 / 12.0      # 4" from O6 east edge to F10
O6_DOOR_WIDTH = 42.0 / 12.0        # 42" door, centered in opening
RO1_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO1
RO2_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO2
RO3_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO3
RO4_DOOR_WIDTH = 36.0 / 12.0       # 36" door in RO4
# Door jamb block thickness = wall - 2*(opening_inside_radius + shell_thickness)
# Opening inside radius 10mm ≈ 0.3937", shell 2", wall 8"
_SHELL = 2.0 / 12.0
_OPENING_R = 10.0 / 304.8
DOOR_FLAT_FACE = WALL_OUTER - 2 * (_OPENING_R + _SHELL)
# O7 (F12-F13 diagonal wall)
O7_NW_GAP = 2.0                    # 2' from F12 to NW end
O7_HALF_WIDTH = 36.0 / 12.0        # 36" half-width (72" total opening)
# O8 (F14-F15)
O8_HALF_WIDTH = 12.5 / 12.0        # 12.5" half-width
# O9, O10, O11
O9_HALF_WIDTH = 12.5 / 12.0        # 12.5" half-width
O10_HALF_WIDTH = 12.5 / 12.0       # 12.5" half-width
O11_HALF_WIDTH = 12.5 / 12.0       # 12.5" half-width

# IW5 partition
IW5_OFFSET_N = 30.0 / 12.0        # 30" south of IW1 south face

# Outline geometry constraints
CORNER_NE_R = 10.0 / 12.0         # R_a0: 10" corner arc
CORNER_NW_R = 28.0 / 12.0         # R_a5: 28" NW corner
UPPER_E_R = 28.0 / 12.0           # R_a7: 28" upper east
SMALL_ARC_R = 2.0 / 12.0          # R_a8: 2" transition
ARC_180_R = 28.0 / 12.0           # R_a11: 28" 180-degree arc
R_a2_a3_DELTA = 8.0 / 12.0        # R_a2 - R_a3 = 8"
F6_HEIGHT = 26.0 - 2.0/12.0       # 25'10" F6-F7 line north of F0
NW_SHIFT = 1.0                    # C5/C3 1' east shift
F1_F2_TARGET = 16.0 + 8.0/12.0    # 16'8" F1-F2 segment target
F4_F5_DROP = 5.0 + 8.0/12.0       # 5'8" F4 south of C5
F16_F17_SEG = 5.0                  # 5' segment
F14_F15_SEG = 8.0 + 4.0/12.0      # 8'4" segment
ARC_F13_R = 5.627004870830987      # R_a13: ~67.52" (set for 75° F10-F11 sweep)
F13_EXIT_BRG = 345.0              # 345-degree exit bearing
SOUTH_WALL_N = -6.0 / 12.0        # -6" south face wall northing
PIX_PI5_TARGET_BRG = 60.0         # 60-degree target bearing
F15_OFFSET_E = 9.0 + 1.0/12.0     # 9'1" F15 east of iw8_e

# Jamb and gap constants
JAMB_WIDTH = 1.0 / 12.0           # 1" jamb width (rough openings)
STD_GAP = 2.0 / 12.0              # 2" standard gap (furniture/appliance spacing)
KITCHEN_APPL_GAP = 3.0 / 12.0     # 3" gap (kitchen appliance spacing/setback)

# Work zone
WW_RADIUS = 30.0 / 12.0           # 30" work-zone radius

# Furniture dimensions
LOVESEAT_WIDTH = 35.0 / 12.0      # 35" loveseat short side
LOVESEAT_LENGTH = 65.0 / 12.0     # 65" loveseat long side
LOVESEAT_ANGLE_DEG = 15.0         # 15° CCW rotation
CHAIR_WIDTH = 32.0 / 12.0         # 32" chair E-W
CHAIR_DEPTH = 37.0 / 12.0         # 37" chair N-S
CHAIR_CORNER_R = 3.0 / 12.0       # 3" rounded corner radius
CHAIR_ANGLE_DEG = 30.0            # 30° CW rotation
OTTOMAN_SIZE = 29.0 / 12.0        # 29" square ottoman
ET_RADIUS_CM = 25.0               # 25 cm endtable radius
SHELVES_WIDTH = 36.0 / 12.0       # 36" shelves E-W
SHELVES_DEPTH = 15.0 / 12.0       # 15" shelves N-S
