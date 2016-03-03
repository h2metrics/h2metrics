function phiPts = evalDiff(pts, phi, splineData)

pts = pts(:);
pts = mod(pts, 2*pi);
phiPts = deBoor( splineData.knotsPhi, splineData.nPhi, phi, pts, ...
                    1, 'periodic', true );
phiPts = phiPts + pts;
phiPts = mod(phiPts, 2*pi);

end