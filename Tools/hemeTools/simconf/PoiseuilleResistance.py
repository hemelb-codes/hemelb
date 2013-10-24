import numpy as np

def SpecificResistance(points, radii):
    """
    Specific => must be multiplied by dynamic viscosity
    
    Assume fully developed Poiseuille flow everywhere.
    
    Define resistance  = pressure drop / flow rate (R = Dp / Q)
    
    Since Poiseuille flow:
    Q = \pi r^4 dp
        ------- --
         8 \eta ds
    where s is the distance along the centreline.
    
    Along each segment assume that the vessel is a conical frustrum with
    
      r_i_____
       |      -------_____r_j
       |                   |
       |                   |
    --s_i-----------------s_j-----> s Centreline
    
    Continuity => Q is a constant.
    So:
      Dp = \int_{s_i}^{s_j} dp ds
                            --
                            ds
    
    Change variables to r(s) = (r_j - r_i) * (s - s_i) + r_i
                               -----------
                               (s_j - s_i)
    
    Dp = 8 \eta Q Ds \int_{r_i}^{r_j} r^{-4} dr
         -------- --
            \pi   Dr
    
       = 8 \eta Q Ds (r_i^{-3} - r_j^{-3})
         -------- --
          3 \pi   Dr
    
    Which is fine unless Dr = r_j - r_i is small
    
    Define: R = (r_i + r_j)/2
    So: r_i = R - Dr/2
    and r_j = R + Dr/2
    
    Expanding in powers of (Dr/(2R))
    r_i^{-3} - r_j^{-3} = 3 Dr   10 Dr^3
                          ---- + -------
                           R^4    4 R^6
    So:
    Dp = 8 \eta Q Ds  /  3    10 Dr^2 \
         ----------- |  --- + -------  |
            3 \pi     \ R^4    4 R^6  /
    
    Let's use this when Dr/R < 0.1
    """
    r_i = radii[:-1]
    r_j = radii[1:]

    Dr = r_j - r_i
    r_mean = 0.5 * (r_i + r_j)
    
    Dx = points[1:] - points[:-1]
    Ds = np.sqrt(np.sum(Dx**2, axis=-1))

    common = (8. * Ds) / (3.0 * np.pi)
    
    simple = common * ( r_i ** -3 - r_j**-3) / Dr
    
    expansion = common * (3.0 * r_mean**-4 +
                          2.5 * Dr**2 * r_mean**-6)

    assert simple.units == expansion.units
    ans = np.where(Dr/r_mean < 0.1, expansion, simple)
    return ans.sum() * simple.units
