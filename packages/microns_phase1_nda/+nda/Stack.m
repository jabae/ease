%{
nda.Stack (manual) # imaging stack
stack_name      : varchar(24)            # name of stack
---
stack_dimensions_x=null     : float                         # 
stack_dimensions_y=null     : float                         # 
stack_dimensions_z=null     : float                         # 
stack_pixels_x=null         : int                           # 
stack_pixels_y=null         : int                           # 
stack_pixels_z=null         : int                           # 
boss_link=null              : varchar(1024)                 # 
stack_description=null      : varchar(1024)                 # free-text description of stack
%}

classdef Stack < dj.Relvar
end