function [K1, K2, phi] = sta_coeffs_distrib(coeffs)
    %STA_COEFFS_DISTRIB Summary of this function goes here
    % build the groups for coeffs of STA
    
    K11 = coeffs(1); K12 = coeffs(2); K13 = coeffs(3);
    K21 = coeffs(4); K22 = coeffs(5); K23 = coeffs(6);
    phi1 = coeffs(7); phi2 = coeffs(8); phi3 = coeffs(9);

    % Start with stronger K1 than K2
    K1 = K12 * ones(9,1);        
    K2 = K22 * ones(9,1); 
    % Smoothing parameter
    phi = phi2 * ones(9,1);
    
    % More aggressive on jet and attitude props
    K1(1) = K11;      
    K2(1) = K21;          
    K1(6:9) = K13;   
    K2(6:9) = K23;      
    
    phi(1) = phi1; 
    phi(6:9) = phi3;  % You may keep this for attitude
end

