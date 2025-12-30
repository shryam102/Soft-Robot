function [u_attr,grad_u_attr] = Attractive(x_goal, x, expo)

u_attr = norm(x - x_goal)^expo;

grad_u_attr = expo*norm(x-x_goal)^(expo - 2)*(x - x_goal);



