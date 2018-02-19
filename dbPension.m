function [ pension ] = dbPension( AIME )
%calibrated to give 25,000 max from 60,000 AIME upwards
%db2 = level/(point^2-2point), db1 = 6.9447e-06*points
%point

    pension =  -(6.9447e-06)*AIME.^2 + 0.8334*AIME;
    pension(AIME > 60000) = 25000;
    
end

