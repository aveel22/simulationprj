import copy
import numpy as np
from typing import Any, Dict


class Parameters(dict):
    """
    A dictionary-like class that allows attribute-style access to parameters.
    Supports both shallow and deep copying with proper handling of nested Parameters objects.
    """

    def __getattr__(self, attr: str) -> Any:
        try:
            return self[attr]
        except KeyError as e:
            raise AttributeError(f"'{type(self).__name__}' has no attribute '{attr}'") from e

    def __setattr__(self, attr: str, value: Any) -> None:
        self[attr] = value

    def __delattr__(self, attr: str) -> None:
        try:
            del self[attr]
        except KeyError as e:
            raise AttributeError(f"'{type(self).__name__}' has no attribute '{attr}'") from e

    def __copy__(self) -> 'Parameters':
        """Create a shallow copy of the Parameters instance."""
        new_obj = Parameters()
        new_obj.update(self)
        return new_obj

    def __deepcopy__(self, memo: Dict[int, Any]) -> 'Parameters':
        """Create a deep copy of the Parameters instance."""
        new_obj = Parameters()
        memo[id(self)] = new_obj  # Add to memo to prevent infinite recursion

        for key, value in self.items():
            # Handle nested Parameters objects specially
            if isinstance(value, Parameters):
                new_obj[key] = copy.deepcopy(value, memo)
            else:
                new_obj[key] = copy.deepcopy(value, memo)

        return new_obj

    def copy(self) -> 'Parameters':
        """Return a shallow copy of the Parameters instance."""
        return copy.copy(self)

    def deepcopy(self) -> 'Parameters':
        """Return a deep copy of the Parameters instance."""
        return copy.deepcopy(self)

    def __variation_geometry(self, scale=1.) -> None:
        nrow, ncol = self.geom_dev.shape
        for i in range(ncol):
            w = self.geom_dev[:, i] * (2 * np.random.rand(nrow) - 1) * scale
            self.geom[:, i] = self.geom[:, i] + w

    def __variate_inertia(self, percent=0.01) -> None:
        nx, ny, nz = self.Inertia.shape
        for i in range(nz):
            self.Inertia[:,:, i] = self.Inertia[:,:, i] * (1. + percent * (2 * np.random.rand(nx, ny) - 1.))

        self.Jo = self.Inertia[:, :, 0]
        self.Jdot *= (1. + percent * (2 * np.random.rand(nx, ny) - 1.))

    def variate_geom(self, scale=1.) -> 'Parameters':
        """
        return copy instance with modified geometry parameters
        """
        new_params:Parameters = self.deepcopy()
        new_params.__variation_geometry(scale)
        return new_params

    def variate_inertia(self, percent=0.01) -> 'Parameters':
        """
        return copy instance with modified inertia parameters
        """
        new_params:Parameters = self.deepcopy()
        new_params.__variate_inertia(percent)
        return new_params

    def variate(self, scale=1., percent=0.) -> None:
        """
        modifies instance parameters
        """
        self.__variate_geom(scale)
        self.__variate_inertia(percent)


def variate_params(param:Parameters, scale=1., percent=0.) -> None:
    """
    creates an instance copy and modifies its parameters
    """
    new_param:Parameters = param.copy()
    return new_param.variate(scale, percent)


def show_matrix(a, name="Matrix A"):
    print(name)
    nx, ny = a.shape
    s = ""
    for i in range(nx):
        for j in range(ny):
            s += f"{a[i, j]:.4}\t"
        s += "\n"
    print(s)


if __name__ == "__main__":
    d = {"k": 4,
         "r": 10,
         "g": 12,
         }
    for k in d.keys():
        print(f'key is {k}, value is {d.get(k)}')
    k = "c"
    print(f'key is {k}, value is {d.get(k)}')
    p2 = Parameters()
    p = Parameters(d)
    print("Parameters:")
    print(p.k)
    print(p["r"])
    variate_params(p, 0.01)