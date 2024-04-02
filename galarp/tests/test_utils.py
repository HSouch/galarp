import numpy as np

import galarp as grp

class TestRotate:
    """ Check that the rotation matrix function actually works as intended, including that it doesn't
        introduce any pesky values when rotating about a principle axis.
    """
    def test_rotate(self):
        # Test rotation of a vector
        xs, ys, zs = np.random.random(100), np.random.random(100), np.random.random(100)
        vec = np.array([xs, ys, zs])
        rotated = grp.rotate(vec, alpha=np.pi/4, beta=np.pi/4, gamma=np.pi/4)
        assert vec.shape == rotated.shape

    def test_rotate_about_axis(self):
        # Test rotation about the z-axis, which should not introduce values in the z-axis (which are all 
        # initialized to 0.)
        
        rand, zero_vec = np.random.random(100), np.zeros(100)

        rotated_x0 = grp.rotate(np.array([zero_vec, rand, rand]), alpha=np.pi/4)
        assert np.max(np.abs(rotated_x0[0])) < 1e-12

        rotated_y0 = grp.rotate(np.array([rand, zero_vec, rand]), beta=np.pi/4)
        assert np.max(np.abs(rotated_y0[1])) < 1e-12

        rotated_z0 = grp.rotate(np.array([rand, rand, zero_vec]), gamma=np.pi/4)
        assert np.max(np.abs(rotated_z0[2])) < 1e-12
