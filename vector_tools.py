from typing import Tuple
import numpy as np 

class BiPoint:
                        #vertices: np.ndarray or Tuple[np.ndarray]
    def __init__(self, vertices: np.ndarray) -> None:
        self.vertices = vertices if isinstance(vertices, np.ndarray) else np.vstack(vertices)
        self._norms = np.linalg.norm(self.vertices, axis=1, keepdims=True) 
        self._normed = self.vertices / self._norms
        self._ip = np.inner(self.vertices[0], self.vertices[1]) 
        self._cp = np.cross(self.vertices[0], self.vertices[1]) 

    def get_midpoint(self) -> np.ndarray:
        return np.sum(self.vertices, axis=0)/2.

    def get_distance(self) -> np.ndarray:
        v = np.diff(self.vertices, axis=0)
        return np.linalg.norm(v)

    def get_angle(self, degree: bool=True) -> float:
        ip = np.inner(self._normed[0], self._normed[1])
        cs = np.clip(ip, -1, 1) # fix rounding issue.
        rads = np.arccos(cs)
        return np.rad2deg(rads) if degree else rads

    def get_cotangent(self) -> float:
        return self._ip/np.linalg.norm(self._cp)
    
    def get_cosecant(self) -> float:
        return np.prod(self._norms)/np.linalg.norm(self._cp)
    
    @classmethod
    def angle(cls, vertices, **kwargs):
        return cls(vertices).get_angle(**kwargs)
    
    @classmethod
    def distance(cls, vertices):
        return cls(vertices).get_distance()
    
    @classmethod
    def midpoint(cls, vertices):
        return cls(vertices).get_midpoint()
    
    @classmethod
    def cotangent(cls, vertices):
        return cls(vertices).get_cotangent()
    
    @classmethod
    def cosecant(cls, vertices):
        return cls(vertices).get_cosecant()
    
class TriPoint:
                        #vertices: np.ndarray or Tuple[np.ndarray]
    def __init__(self, vertices: np.ndarray ) -> None:
        self.vertices = vertices if isinstance(vertices, np.ndarray) else np.vstack(vertices)
        self._lengths = self.get_lenghts()
        self._vectors = self.vertices[1:] - self.vertices[0][np.newaxis, :]

    def get_area(self) -> float:
        # https://en.wikipedia.org/wiki/Heron%27s_formula
        s = np.sum(self._lengths)/2.
        mask = np.isclose(s, self._lengths) # fix rounding issue.
        aux = np.where(mask, 0., s-self._lengths)
        return np.sqrt(s*np.prod(aux))

    def get_lenghts(self) -> np.ndarray:
        vectors = np.diff(self.vertices, axis=0, append=self.vertices[0][np.newaxis,:])
        return np.linalg.norm(vectors, axis=1)

    def get_compactness(self) -> float:
        return 4*np.sqrt(3.)*self.get_area()/np.linalg.norm(self._lengths)**2

    def get_normal(self, normalize: bool = True) -> np.ndarray:
        n = np.cross(self._vectors[0], self._vectors[1])
        if normalize: n /= np.linalg.norm(n)
        return n

    def get_angle(self, **kwargs) -> float:
        return BiPoint.angle(self._vectors, **kwargs) 
    
    def get_cotangent(self) -> float:
        return BiPoint.cotangent(self._vectors)
    
    def get_cosecant(self) -> float:
        return BiPoint.cosecant(self._vectors)
    
    def get_is_collinear(self) -> bool:
        n = self.get_normal(normalize=False)
        return np.allclose(n, 0)
    
    def get_is_obtuse(self) -> bool:
        # https://mathworld.wolfram.com/LawofCosines.html
        # https://mathworld.wolfram.com/ObtuseTriangle.html
        a, b, c = self._lengths[0], self._lengths[1], self._lengths[2] 
        return (a**2+c**2 < b**2) or (a**2+b**2 < c**2) or (b**2+c**2 < a**2)
    
    def get_barycentric_region(self) -> float:
        return self.get_area()/3. 

    def get_voronoi_region(self) -> float:
        a, _, c = self._lengths[0], self._lengths[1], self._lengths[2] 
        p, r, q = self.vertices[0], self.vertices[1], self.vertices[2] 
        return (1/8.) * (a**2 * BiPoint.cotangent((p-q, r-q)) + c**2 * BiPoint.cotangent((p-r, q-r)))
    
    def get_mixed_voronoi_region(self) -> float:
        # http://rodolphe-vaillant.fr/entry/20/compute-harmonic-weights-on-a-triangular-mesh#mixed_voro_area
        if not self.get_is_obtuse(): return self.get_voronoi_region()
        a, b, c = self._lengths[0], self._lengths[1], self._lengths[2]
        return self.get_area()/2. if (a**2+c**2 < b**2) else self.get_area()/4.

    def get_barycenter(self) -> np.ndarray:
        return np.mean(self.vertices, axis=0)
    
    def get_barycentric_coords(self, point: np.ndarray) -> np.ndarray:
        p, r, q = self.vertices[0], self.vertices[1], self.vertices[2] 
        coords = np.array([TriPoint.area((point, q, r)), TriPoint.area((point, p, q)), TriPoint.area((point, p, r))])
        coords /= np.sum(coords)
        return coords
    
    def get_barycentric_interp(self, values: np.ndarray) -> float:
        bary = self.get_barycenter()
        coords = self.get_barycentric_coords(bary)
        return np.inner(values, coords)
    
    @classmethod
    def is_collinear(cls, vertices):
        return cls(vertices).get_is_collinear()

    @classmethod
    def barycenter(cls, vertices):
        return cls(vertices).get_barycenter()
    
    @classmethod
    def normal(cls, vertices, **kwargs):
        return cls(vertices).get_normal(**kwargs)
    
    @classmethod
    def area(cls, vertices):
        return cls(vertices).get_area()
    
    @classmethod
    def barycentric_interp(cls, vertices, *args, **kwargs):
        return cls(vertices).get_barycentric_interp(*args, **kwargs)