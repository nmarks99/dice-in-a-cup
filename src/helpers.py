import sympy as sym
import numpy as np
from sympy.core import symbol


def SE3(p,theta,symbolic=True):
    if symbolic:
        g = sym.Matrix([
            [sym.cos(theta), -sym.sin(theta),    0, p[0]],
            [sym.sin(theta),  sym.cos(theta),    0, p[1]],
            [             0,               0,  1.0,    0],
            [             0,               0,    0,  1.0]
        ])
    
    else:
        g = np.array([
            [np.cos(theta), -np.sin(theta),    0, p[0]],
            [np.sin(theta),  np.cos(theta),    0, p[1]],
            [            0,              0,  1.0,    0],
            [            0,              0,    0,  1.0]
        ])
    return g 



def hat(w,symbolic=True):
    ''' 
    hat() does the "hat" operation which generates a 3x3 skew-symmetric matrix given a 1x3 vector 
    '''

    if symbolic:
        w_hat = sym.Matrix(
            [
                [    0, -w[2],  w[1]],
                [ w[2],     0, -w[0]],
                [-w[1],  w[0],     0]
            ])
    else:
        w_hat = np.array(
            [
                [    0, -w[2],  w[1]],
                [ w[2],     0, -w[0]],
                [-w[1],  w[0],     0]
            ])
    return w_hat


def unhat(mat):
    '''
    unhat() does the unhat operation on a 4x4 matrix which is in our case often just the product of 
    g_inv and g_dot for use in calculating the body velocity
    '''
    xdim,ydim = mat.shape
    assert(xdim == ydim),'Input must be square 4x4 matrix'
    assert(xdim == 4),'Input must be 4x4 matrix'
    vec = sym.Matrix([mat[0,3],mat[1,3],mat[2,3],mat[2,1],mat[0,2],mat[1,0]])
    return vec
    

def invert_g(g):
    '''
    invert_g takes finds the inverse of the transformation matrix g by indexing instead of 
    having to actually compute the inverse using a standard method. This is faster.
    '''
    assert(g.shape == (4,4)),'This is only applicable to 4x4 transformation matrices'

    zeros_3x1 = sym.Matrix([0,0,0])
    one = sym.Matrix([1])
    
    R = sym.Matrix(
        [
            [g[0,0],g[0,1],g[0,2]], 
            [g[1,0],g[1,1],g[1,2]], 
            [g[2,0],g[2,1],g[2,2]]
        ])
    p = sym.Matrix(
        [
            g[0,3],g[1,3],g[2,3]
        ])
    
    g_inv = sym.Matrix(
        [
            [R.T, -R.T*p],
            [zeros_3x1.T, one]            
        ])
        
    return g_inv

def get_Vb(g):
    '''
    get_Vb(g) computes the body velocity of the inputted tranformation matrix by doing unhat(g_inv * g_dot)
    '''
    from sympy.abc import t
    g_inv = invert_g(g)
    g_dot = g.diff(t)
    Vb = unhat(g_inv*g_dot)
    return Vb