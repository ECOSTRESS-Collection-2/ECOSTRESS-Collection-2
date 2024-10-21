class Tower:
    def __init__(self, name: str, lat: float, lon: float):
        self.name = name
        self.lat = lat
        self.lon = lon

    def __repr__(self) -> str:
        return f'Tower(name="{self.name}", lat={self.lat}, lon={self.lon})'
