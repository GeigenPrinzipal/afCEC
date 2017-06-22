afCEC <- function(points, maxClusters, initialLabels, cardMin, minIterations, maxIterations, numberOfStarts, method, values, interactive=FALSE) {
    if (class(values) == "character") {
        if (values == "quadratic") {
            pointsNum = dim(points)[1];
            dimension = dim(points)[2];
            valuesArray = matrix(rep(0,((2*(dimension-1))+1)*dimension*pointsNum),((2*(dimension-1))+1)*dimension,pointsNum);
            for (i in 1:pointsNum) {
                tmp = points[i,2:dimension];
                for (j in 1:dimension) {
                    for (k in 1:(dimension - 1)) {
                        valuesArray[((j-1)*((2*(dimension-1))+1))+k,i]=tmp[k]^2;
                        valuesArray[((j-1)*((2*(dimension-1))+1))+dimension-1+k,i]=tmp[k];
                        valuesArray[((j-1)*((2*(dimension-1))+1))+(2*(dimension-1))+1,i]=1;
                    }
                    if (j < dimension) tmp[j] = points[i, j];
                }
            }
        } else {
            CalculateArrayOfValuesCode = GenerateCodeForArrayConstruction(values);
            writeLines(CalculateArrayOfValuesCode);
            sourceCpp(code=CalculateArrayOfValuesCode);
            valuesArray = CalculateArrayOfValues(t(points));
        }
    } else {
        valuesArray = values;
    }
    res=afCECCppRoutine(t(points),maxClusters,initialLabels,cardMin,minIterations,maxIterations,numberOfStarts,method,valuesArray,interactive);
    res[["values"]]=valuesArray;
    if (length(res) > 0) {
        if (class(values) == "character") {
            if (values != "quadratic") {
                UpdateMeansCode = GenerateCodeForUpdatingMeans(values);
                writeLines(UpdateMeansCode);
                sourceCpp(code=UpdateMeansCode);
                if (!interactive) {
                    UpdateMeans(res);
                    res[["data"]] = points;
                    res[["formula"]] = values;
                } else {
                    for (i in 1:length(res)) {
                        UpdateMeans(res[[i]]);
                        res[[i]][["data"]] = points;
                        res[[i]][["formula"]] = values;
                    }
                }
            } else {
                if (!interactive) {
                    UpdateMeansForQuadraticFunction(res);
                    res[["data"]] = points;
                    res[["formula"]] = values;
                } else {
                    for (i in 1:length(res)) {
                        UpdateMeansForQuadraticFunction(res[[i]]);
                        res[[i]][["data"]] = points;
                        res[[i]][["formula"]] = values;
                    }
                }
            }
        } else {
            if (!interactive) {
                res[["data"]] = points;
                res[["formula"]] = NULL;
            } else {
                for (i in 1:length(res)) {
                    res[[i]][["data"]] = points;
                    res[[i]][["formula"]] = NULL;
                }
            }
        }
    }
    if (!interactive) {
        res=structure(res,class="afCEC");
    } else {
        for (i in 1:length(res)) res[[i]] = structure(res[[i]], class="afCEC");
    }
    return(res);
}

#   -*-   -*-   -*-

plot.afCEC <- function(x, draw_points=TRUE, draw_means=TRUE, draw_ellipsoids=TRUE, draw_surfaces=TRUE, confidence=0.95, grid_resolution=32) {
    if (dim(x$data)[2] == 2) {
        if (draw_points) plot(x$data,pch=20,col=x$labels+1,asp=1)
        if ((!is.null(x$formula)) && (draw_means || draw_ellipsoids || draw_surfaces)) {
            if (draw_ellipsoids || draw_surfaces) {
                if (x$formula != "quadratic") {
                    # formula
                    CalculateEllipsesOfConfidenceCode = GenerateCodeForCalculatingEllipseOfConfidence(x$formula);
                    sourceCpp(code=CalculateEllipsesOfConfidenceCode);
                    ellipses=CalculateEllipsesOfConfidence(x, confidence, grid_resolution);
                } else {
                    # quadratic function
                    ellipses=CalculateEllipsesOfConfidenceForQuadraticFunction(x, confidence, grid_resolution);
                }
            }
            for (i in 1:length(x$means)) {
                if (!is.null(x$means[[i]])) {
                    if (draw_means) points(t(x$means[[i]]),pch=20,cex=2.5,col=rgb(0,0,0));
                    if (draw_ellipsoids) lines(ellipses[[1]][[i]]);
                    if (draw_surfaces) lines(ellipses[[2]][[i]]);
                }
            }
        }
    } else {
        if (dim(x$data)[2] == 3) {
            if (draw_points) plot3d(x$data,col=x$labels+1,aspect=FALSE,alpha=1,xlab="",ylab="",zlab="");
            if ((!is.null(x$formula)) && (draw_means || draw_ellipsoids || draw_surfaces)) {
                if (draw_ellipsoids || draw_surfaces) {
                    if (x$formula != "quadratic") {
                        #formula
                        CalculateEllipsoidsOfConfidenceCode = GenerateCodeForCalculatingEllipsoidsOfConfidence(x$formula);
                        sourceCpp(code=CalculateEllipsoidsOfConfidenceCode);
                        ellipsoids=CalculateEllipsoidsOfConfidence(x, confidence, grid_resolution);
                    } else {
                        #quadratic function (to be implemented for 3d data)
                    }
                }
                if (draw_means) {
                    dataSize=dim(x$data)[1];
                    dimX=max(x$data[1:dataSize,1])-min(x$data[1:dataSize,1]);
                    dimY=max(x$data[1:dataSize,2])-min(x$data[1:dataSize,2]);
                    dimZ=max(x$data[1:dataSize,3])-min(x$data[1:dataSize,3]);
                    maxDim=max(dimX,dimY,dimZ);
                }
                for (i in 1:length(x$means)) {
                    if (!is.null(x$means[[i]])) {
                        if (draw_ellipsoids) triangles3d(ellipsoids[[1]][[i]],normals=ellipsoids[[2]][[i]],alpha=0.25,col=i);
                        if (draw_surfaces) triangles3d(ellipsoids[[3]][[i]],normals=ellipsoids[[4]][[i]],alpha=0.5,col=i);
                        if (draw_means) spheres3d(t(x$means[[i]]),radius=maxDim*0.01,col=rgb(0,0,0),alpha=0.5);
                    }
                }
            }
        }
    }
}

#   -*-   -*-   -*-

#####################################
#                                   #
# Functions for data stream testing #
#                                   #
#####################################

addPoint <- function(x, point) UseMethod("addPoint");

# x - object of the type afCEC
# point - vector with a coordinates of the point to be added
addPoint.afCEC <- function(x, point) {
    n = dim(x$data)[1];
    d = dim(x$data)[2];
    cl = length(x$cardinalities);
    # macierz sumy kwadratów różnych współrzędnych punktów należących do poszczególnych klastrów pozwoli na sprawne
    # wyznaczenie wartości wariancji na osi aktywnej
    sumOfSquares = matrix(rep(0,cl*d),cl,d);

    for (i in 1:cl) {
        for (j in 1:d) {
            sumOfSquares[i,j]=0;
        }
    }
    for (i in 1:cl) {
        if (x$cardinalities[[i]] > 0) {
            x$covariances[[i]] = matrix(rep(0,d*d),d,d);
            x$means[[i]] = matrix(rep(0,d),d,1);
        } else {
            x$covariances[[i]] = NULL;
            x$means[[i]] = NULL;
        }
    }

    for (i in 1:n) {
        lbl = x$labels[i] + 1;
        x$means[[lbl]] = x$means[[lbl]]+x$data[i,];
    }
    for (i in 1:cl) {
        if (x$cardinalities[[i]] > 0) {
            x$means[[i]] = x$means[[i]]/x$cardinalities[[i]];
        }
    }

    for (i in 1:n) {
        lbl = x$labels[i] + 1;
        for (j in 1:d) {
            sumOfSquares[lbl,j]=sumOfSquares[lbl,j]+x$data[i,j]^2;
        }
        x$covariances[[lbl]] = x$covariances[[lbl]]+((x$data[i,]-x$means[[lbl]])%*%t(x$data[i,]-x$means[[lbl]]));
    }
    for (i in 1:cl) {
        if (x$cardinalities[[i]] > 0) {
            x$covariances[[i]] = x$covariances[[i]]/x$cardinalities[[i]];
        }
    }

    dCostMin = Inf;
    EMin = NA;
    for (i in 1:cl) {
        if (x$cardinalities[[i]] > 0) {
            # wektor tmp w każdym kroku pętli po j przechowuje współrzędne punktu [point] z wyrzuconą [j]-tą współrzędną
            tmp=point[2:d];
            for (j in 1:d) {
                tryCatch({
                    covTmp=(x$cardinalities[[i]]/(x$cardinalities[[i]]+1))*
                        (x$covariances[[i]]+((1/(x$cardinalities[[i]]+1))*((point-x$means[[i]])%*%t(point-x$means[[i]]))));

                    # Wyznaczamy wartość funkcji aktywnej w punkcie [point] dla osi aktywnej [j]
                    val=rep(0,(2*(d-1))+1);
                    for (k in 1:(d-1)) {
                        val[k]=tmp[k]^2;
                        val[d-1+k]=tmp[k];
                    }
                    val[(2*(d-1))+1]=1;
                    if (j < d) tmp[j]=point[j];

                    ATmp=x$GLM_design_matrix[[i]][[j]]+(val%*%t(val));
                    bTmp=x$GLM_response_vector[[i]][[j]]+(point[j]*val);
                    coeffs=solve(ATmp,bTmp);
                    var=sumOfSquares[i,j]+(point[j]^2);
                    var=var-(2*(t(coeffs)%*%bTmp));
                    var=var+(t(coeffs)%*%ATmp%*%coeffs);
                    var=var/(x$cardinalities[[i]]+1);

                    prob=(x$cardinalities[[i]]+1)/n;
                    detVal = det(as.matrix(covTmp[-j,-j]))*var;
                    if (detVal < 0.0) detVal = 0.0;
                    ETmp=prob*(-log(prob)+(0.5*((d*log(2*pi*exp(1)))+log(detVal))));
                    if ((is.finite(ETmp)) && (is.finite(x$cost[[i]]))) {
                        if (ETmp-x$cost[[i]] < dCostMin) {
                            dCostMin=ETmp-x$cost[[i]];
                            bestCl=i;
                            # !!!
                            bestDir=j;
                            bestCov=covTmp;
                            bestVar=var;
                            bestCoeffs=coeffs;
                            # !!!
                            EMin=ETmp;
                        }
                    }
                }, error=function(e) {
                });
            }
        }
    }
    tmp=point[2:d];
    for (i in 1:d) {
        val=rep(0,(2*(d-1))+1);
        for (j in 1:(d-1)) {
            val[j]=tmp[j]^2;
            val[d-1+j]=tmp[j];
        }
        val[(2*(d-1))+1]=1;
        if (i < d) tmp[i]=point[i];

        x$GLM_design_matrix[[bestCl]][[i]]=x$GLM_design_matrix[[bestCl]][[i]]+(val%*%t(val));
        x$GLM_response_vector[[bestCl]][[i]]=x$GLM_response_vector[[bestCl]][[i]]+(point[i]*val);
    }

    x$means[[bestCl]]=((x$means[[bestCl]]*x$cardinalities[[bestCl]])+point)/(x$cardinalities[[bestCl]]+1);

    # !!!
    x$covariances[[bestCl]]=bestCov;
    if (bestDir>1) {
        x$covariances[[bestCl]][bestDir,1:(bestDir-1)]=0;
        x$covariances[[bestCl]][1:(bestDir-1),bestDir]=0;
    }
    if (bestDir<d) {
        x$covariances[[bestCl]][bestDir,(bestDir+1):d]=0;
        x$covariances[[bestCl]][(bestDir+1):d,bestDir]=0;
    }
    x$covariances[[bestCl]][bestDir,bestDir]=bestVar;
    x$directions[[bestCl]]=bestDir-1;
    x$coefficients[[bestCl]]=bestCoeffs;
    # !!!

    x$cardinalities[[bestCl]]=x$cardinalities[[bestCl]]+1;
    x$cost[[bestCl]]=EMin;
    x$data=rbind(x$data,point);
    x$labels=c(x$labels,bestCl-1);
    return(x);
}

#   -*-   -*-   -*-

removePoint <- function(x, point) UseMethod("removePoint");

# x - object of the type afCEC
# ind - index of the point to be removed
removePoint.afCEC <- function(x, ind) {
    n = dim(x$data)[1];
    d = dim(x$data)[2];
    cl = length(x$cardinalities);
    # macierz sumy kwadratów różnych współrzędnych punktów należących do poszczególnych klastrów pozwoli na sprawne
    # wyznaczenie wartości wariancji na osi aktywnej
    sumOfSquares = matrix(rep(0,cl*d),cl,d);

    for (i in 1:cl) {
        for (j in 1:d) {
            sumOfSquares[i,j]=0;
        }
    }
    for (i in 1:cl) {
        if (x$cardinalities[[i]] > 0) {
            x$covariances[[i]] = matrix(rep(0,d*d),d,d);
            x$means[[i]] = matrix(rep(0,d),d,1);
        } else {
            x$covariances[[i]] = NULL;
            x$means[[i]] = NULL;
        }
    }

    for (i in 1:n) {
        lbl = x$labels[i] + 1;
        x$means[[lbl]] = x$means[[lbl]]+x$data[i,];
    }
    for (i in 1:cl) {
        if (x$cardinalities[[i]] > 0) {
            x$means[[i]] = x$means[[i]]/x$cardinalities[[i]];
        }
    }

    for (i in 1:n) {
        lbl = x$labels[i] + 1;
        for (j in 1:d) {
            sumOfSquares[lbl,j]=sumOfSquares[lbl,j]+x$data[i,j]^2;
        }
        x$covariances[[lbl]] = x$covariances[[lbl]]+((x$data[i,]-x$means[[lbl]])%*%t(x$data[i,]-x$means[[lbl]]));
    }
    for (i in 1:cl) {
        if (x$cardinalities[[i]] > 0) {
            x$covariances[[i]] = x$covariances[[i]]/x$cardinalities[[i]];
        }
    }

    # usuwany punkt
    point=x$data[ind,];
    # etykieta usuwanego punktu
    lbl=x$labels[ind]+1;
    # wektor tmp w każdym kroku pętli po i przechowuje współrzędne punktu [point] z wyrzuconą [i]-tą współrzędną
    tmp=point[2:d];
    dCostMin = Inf;
    EMin = NA;
    for (i in 1:d) {
        tryCatch({
            covTmp=(x$cardinalities[[lbl]]/(x$cardinalities[[lbl]]-1))*
                (x$covariances[[lbl]]-((1/(x$cardinalities[[lbl]]-1))*((point-x$means[[lbl]])%*%t(point-x$means[[lbl]]))));

            # Wyznaczamy wartość funkcji aktywnej w punkcie [point] dla osi aktywnej [i]
            val=rep(0,(2*(d-1))+1);
            for (j in 1:(d-1)) {
                val[j]=tmp[j]^2;
                val[d-1+j]=tmp[j];
            }
            val[(2*(d-1))+1]=1;
            if (i < d) tmp[i]=point[i];

            ATmp=x$GLM_design_matrix[[lbl]][[i]]-(val%*%t(val));
            bTmp=x$GLM_response_vector[[lbl]][[i]]-(point[i]*val);
            coeffs=solve(ATmp,bTmp);
            var=sumOfSquares[lbl,i]-(point[i]^2);
            var=var-(2*(t(coeffs)%*%bTmp));
            var=var+(t(coeffs)%*%ATmp%*%coeffs);
            var=var/(x$cardinalities[[lbl]]-1);

            prob=(x$cardinalities[[lbl]]-1)/n;
            detVal = det(as.matrix(covTmp[-i,-i]))*var;
            if (detVal < 0.0) detVal = 0.0;
            ETmp=prob*(-log(prob)+(0.5*((d*log(2*pi*exp(1)))+log(detVal))));
            if ((is.finite(ETmp)) && (is.finite(x$cost[[lbl]]))) {
                if (ETmp-x$cost[[lbl]] < dCostMin) {
                    dCostMin=ETmp-x$cost[[lbl]];
                    # !!!
                    bestDir=i;
                    bestCov=covTmp;
                    bestVar=var;
                    bestCoeffs=coeffs;
                    # !!!
                    EMin=ETmp;
                }
            }
        }, error=function(e) {
        });
    }
    tmp=point[2:d];
    for (i in 1:d) {
        val=rep(0,(2*(d-1))+1);
        for (j in 1:(d-1)) {
            val[j]=tmp[j]^2;
            val[d-1+j]=tmp[j];
        }
        val[(2*(d-1))+1]=1;
        if (i < d) tmp[i]=point[i];

        x$GLM_design_matrix[[lbl]][[i]]=x$GLM_design_matrix[[lbl]][[i]]-(val%*%t(val));
        x$GLM_response_vector[[lbl]][[i]]=x$GLM_response_vector[[lbl]][[i]]-(point[i]*val);
    }
    x$means[[lbl]]=((x$means[[lbl]]*x$cardinalities[[lbl]])-point)/(x$cardinalities[[lbl]]-1);

    # !!!
    if (!is.na(EMin)) {
        x$covariances[[lbl]]=bestCov;
        if (bestDir>1) {
            x$covariances[[lbl]][bestDir,1:(bestDir-1)]=0;
            x$covariances[[lbl]][1:(bestDir-1),bestDir]=0;
        }
        if (bestDir<d) {
            x$covariances[[lbl]][bestDir,(bestDir+1):d]=0;
            x$covariances[[lbl]][(bestDir+1):d,bestDir]=0;
        }
        x$covariances[[lbl]][bestDir,bestDir]=bestVar;
        x$directions[[lbl]]=bestDir-1;
        x$coefficients[[lbl]]=bestCoeffs;
    }
    # !!!

    x$cardinalities[[lbl]]=x$cardinalities[[lbl]]-1;
    x$cost[[lbl]]=EMin;
    x$data=x$data[-ind,];
    x$labels=x$labels[-ind];
    return(x);
}

#   -*-   -*-   -*-

SampleTorusUniform <- function(R, r, numberOfPoints) {
    points=matrix(rep(0, numberOfPoints*3), numberOfPoints, 3);
    for (i in 1:numberOfPoints) {
        theta=2*pi*runif(1, 0, 1);
        phi=2*pi*runif(1, 0, 1);
        x=(R+(r*cos(phi)))*cos(theta);
        y=(R+(r*cos(phi)))*sin(theta);
        z=r*sin(phi);
        points[i, 1]=x;
        points[i, 2]=y;
        points[i, 3]=z;
    }
    return(points);
}
