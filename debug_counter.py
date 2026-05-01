import json
from app.engine import compute_geometry
from app.database import get_db, get_constants_dict, get_outline_chain, get_all_doors

variant = 'daybed'
constants = get_constants_dict('app/adu.db')
chain = get_outline_chain('app/adu.db')
doors = get_all_doors('app/adu.db')
geom = compute_geometry(constants, variant, chain, doors_data=doors, db_path='app/adu.db')
print('variant_items count', len(geom.get('variant_items', {})))
for name in ['COUNTER', 'counter', 'STOVE', 'MICRO', 'microwave', 'kitchen_counter']:
    if name in geom.get('variant_items', {}):
        item = geom['variant_items'][name]
        print('FOUND', name, 'type', item.get('type'), 'shape', item.get('shape'), 'poly_len', len(item.get('poly', [])), 'center', item.get('center'))

with get_db() as conn:
    rows = conn.execute("SELECT id,name,type,variant,properties FROM elements WHERE type IN ('furniture','appliance','fixture') ORDER BY name").fetchall()
    print('DB items count', len(rows))
    for r in rows:
        name = r['name']
        if name.lower() in ('counter', 'stove', 'micro', 'kitchen_counter', 'custom_daybed') or name in ('COUNTER', 'STOVE', 'MICRO'):
            props = r['properties']
            print('DB', dict(r), 'props_center', 'center' in props)
    print('All matching DB names:')
    for r in rows:
        name = r['name']
        if any(k in name.lower() for k in ('counter', 'stove', 'micro')):
            print('  ', dict(r))
