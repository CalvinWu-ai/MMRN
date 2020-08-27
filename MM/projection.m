function Up = projection(X, U)
        
            
% The projection onto the tangent space of the Stiefel manifold.
%         
    Up = U - X*symm(X'*U);
   
 end