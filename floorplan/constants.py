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
KITCHEN_SINK_WIDTH = 48.0 / 12.0   # 48" E-W
KITCHEN_SINK_DEPTH = 24.0 / 12.0   # 24" N-S
DW_WIDTH = 28.0 / 12.0             # 28" dishwasher E-W
DW_DEPTH = 27.0 / 12.0             # 27" dishwasher N-S
STOVE_WIDTH = 30.0 / 12.0          # 30" stove E-W
STOVE_DEPTH = 27.0 / 12.0          # 27" stove N-S
FRIDGE_SIZE = 36.0 / 12.0          # 36" fridge (square)
KITCHEN_GAP = 0.75 / 12.0          # 3/4" gap between kitchen appliances
KITCHEN_CTR_LENGTH = 72.0 / 12.0   # 72" kitchen counter E-W along IW1 north
KITCHEN_CTR_DEPTH = 24.0 / 12.0    # 24" kitchen counter depth N-S
NORTH_CTR_LENGTH = 38.0 / 12.0     # 38" north wall counter E-W
NORTH_CTR_DEPTH = 24.0 / 12.0      # 24" north wall counter depth N-S

# IW6 partition
IW6_THICKNESS = 1.0 / 12.0        # 1" partition
IW6_OFFSET_N = 5.5                 # 5'6" south of F6-F7 south face

# IW5 partition
IW5_OFFSET_N = 30.0 / 12.0        # 30" south of IW1 south face

# Outline geometry constraints
CORNER_NE_R = 10.0 / 12.0         # R_a0: 10" corner arc
CORNER_NW_R = 28.0 / 12.0         # R_a5: 28" NW corner
UPPER_E_R = 28.0 / 12.0           # R_a7: 28" upper east
SMALL_ARC_R = 2.0 / 12.0          # R_a8, R_a10: 2" transitions
ARC_180_R = 28.0 / 12.0           # R_a11: 28" 180-degree arc
R_a2_a3_DELTA = 8.0 / 12.0        # R_a2 - R_a3 = 8"
F6_HEIGHT = 26.0 - 2.0/12.0       # 25'10" F6-F7 line north of F0
NW_SHIFT = 1.0                    # C5/C3 1' east shift
F1_F2_TARGET = 16.0 + 8.0/12.0    # 16'8" F1-F2 segment target
F4_F5_DROP = 5.0 + 8.0/12.0       # 5'8" F4 south of C5
F16_F17_SEG = 5.0                  # 5' segment
F14_F15_SEG = 8.0 + 4.0/12.0      # 8'4" segment
F12_F13_SEG = 10.0                 # 10' segment
F13_EXIT_BRG = 345.0              # 345-degree exit bearing
SOUTH_WALL_N = -6.0 / 12.0        # -6" south face wall northing
PIX_PI5_TARGET_BRG = 60.0         # 60-degree target bearing
F15_OFFSET_E = 9.0 + 1.0/12.0     # 9'1" F15 east of iw8_e
