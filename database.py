import json
def load_solar_system(filepath='solar_system.json'):
    with open(filepath,'r',encoding='utf-8') as f:
        data = json.load(f)
        return data