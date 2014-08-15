"""
Stores all the colors for the main theme.
"""

## Blue colors
#bc = '#006C98'
bc = '#008FC0'
lbc = '#3292CB'#  		'#3B9BC4'
dbc = '#00364C'

bcscr = '#940014'
bccr = '#942A00'		# Reddish color

## Red colors
#rc = '#D20202'
rc = '#CF363A'
drc = '#690101'
lrc = '#E88080'			# Pink

rcscb = '#027ACF'

## Green colors
#gc = '#007935';	
gc = '#40993B';									
lgc = '#69f432';
dgc = '#145114';

# Green complements
gcscp = '#702AA2';			# Split complement - purple
gcscb = '#2A48A2';			# Split complement - blue
gctcr = '#A22A2A';			# Triadic color - red
gctcb = '#2A2AA2';			# Triadic color - blue
gccp = '#A22AA2';			# Complement - purple

## Purple colors
#pc = '#AB4CEE'
pc = '#84439A'
lpc = '#E3A0EE'
dpc = '#552677'

## Cyan colors
cc = '#3AA9AA'
lcc = '#8CD7D8'
dcc = '#1A4B4C'

## Orange colors
oc = '#ffa71a'
loc = 'ffc05d'				# Tan
doc = 'd44f00'

## Grey colors
grc = '#A29F96'
lgrc = '#D1CFCB'
dgrc = '#3C3B38'

## Yellow colors
yc = '#fdef19'
lyc = '#FFFF66'
dyc = '#A37B18'


# Color cycles
color_cycle = [bc, rc, gc, oc, grc, pc, yc, cc, dgrc]
dark_color_cycle = [dbc, drc, dgc, dpc, dyc, dcc, doc]
light_color_cycle = [lbc, lrc, lgc, lpc, lcc, loc, lgrc, lyc]

extended_color_cycle = [bc, rc, gc, oc, grc, pc, yc, cc, dgrc, gcscp, gctcr, gctcb, lrc, lbc, lgc, lpc, lcc, loc, lgrc, dbc, drc, dgc, dpc, dyc, dcc, doc]

green_complements_color_cycle = [gc, gccp, gctcr, gctcb, grc, gcscp, gcscb]