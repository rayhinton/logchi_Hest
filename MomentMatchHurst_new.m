function [h_hat, E] = MomentMatchHurst_new(wddata, pairs, L, ismean)
    
    J = log2(size(wddata,2));

    H_hat = zeros(size(pairs, 1), 2);
    
    E = [];
    for p = 1: size(pairs,1)
                
        k1 = pairs(p,1); k2 = pairs(p,2); 

        %l1 =  J - k1; l2 =  J - k2;
        % l1 =  round(k1-L+1); l2 =  round(k2-L+1);
        % RJH k1 and k2 do not need to be converted; they can be used as-is
        % to get calculate the indices for wddata.
        l1 = k1; l2 = k2;
        k1_indx =  2^(l1) + 1 : 2^(l1 + 1); 
        k2_indx =  2^(l2) + 1 : 2^(l2 + 1);
                
        %length(k1_indx) length(k2_indx)] 

        d_k1 = log2( median ( wddata( k1_indx ).^2));  
        d_k2 = log2( median ( wddata( k2_indx ).^2));

        A = psi( 2^( l2 - 1 ) ) - psi( 2^(l1 - 1) ) ;
        B = log2( exp(1) ) * A - d_k2 + d_k1;
        
        C = l2 - l1;
        h_p = ( 1/ ( 2*C ) ) * ( B  ) - 1; 

        h_p_var = ( 1/ ( (2*C)^2 ) ) * ( psi( 1, 2^( l2 - 2 ) ) - psi( 1, 2^(l1 - 2) ) );
                
        H_hat(p, :) =  [ h_p h_p_var ];
        E = [E; [d_k1 d_k2]];
   end
            
   wgt = H_hat(:, 2)/sum(H_hat(:, 2));
   if ismean == 1        
    h_hat = mean(H_hat(:,1));% 
    %h_hat =  H_hat(:,1)'* wgt;
   elseif ismean == 0
    h_hat = median(H_hat(:,1));
    %h_hat = H_hat(:,1)'* wgt;
   end