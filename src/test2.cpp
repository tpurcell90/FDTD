if(Ez_)
    {
        for(int kk = 0; kk < pmlArr_.size(); kk++)
        {
            switch(pmlArr_[kk].d())
            {
                case X:
                {
                    if(yPML_ ==0)
                        a = 0;
                    else
                        a = yPML_;
                    for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                    {
                        double eps = 1.0;
                        double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                        double sigz = 0.0;
                        double sigxx = pmlArr_[kk].sigma(static_cast<double>(ii));
                        double sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5);
                        double sigyx = 0.0;
                        double sigyy = 0.0;
                       
                        double c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                        double c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;
                        double c_hxh = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                        double c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                        double c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                        double c_byb = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                        double c_bye = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_)) ;
                        double c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                        double c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                        double c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                     
                        for(int jj =  a; jj < ny_ - a; jj ++)
                        {

                            double bxstore = pmlArr_[kk].Bx_->point(ii,jj);
                            double bystore = pmlArr_[kk].By_->point(ii,jj);
                            if(jj != ny_-1)
                            {
                                pmlArr_[kk].Bx_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                pmlArr_[kk].By_->point(ii,jj) = c_byb * pmlArr_[kk].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                                
                                Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                                Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[kk].By_->point(ii,jj) - c_hyb0 * bystore;
                                
                                sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) - 0.5);
                                c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                
                                bxstore = pmlArr_[kk].Bx_end_->point(ii,jj);
                                bystore = pmlArr_[kk].By_end_->point(ii,jj);

                                pmlArr_[kk].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point((nx_-1) - ii, jj+1)-Ez_->point((nx_-1) - ii,jj));
                                pmlArr_[kk].By_end_->point(ii,jj) = c_byb * pmlArr_[kk].By_end_->point(ii,jj) + c_bye * (Ez_->point((nx_-1) - ii+1, jj)-Ez_->point((nx_-1) - ii,jj));

                                Hx_->point((nx_-1) - ii,jj) = c_hxh * Hx_->point((nx_-1) - ii,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                                Hy_->point((nx_-1) - ii,jj) = c_hyh * Hy_->point((nx_-1) - ii,jj) + c_hyb1 * pmlArr_[kk].By_end_->point(ii,jj) - c_hyb0 * bystore;
                            }
                            else
                            {
                                pmlArr_[kk].By_->point(ii,jj) = c_byb * pmlArr_[kk].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                                Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[kk].By_->point(ii,jj) - c_hyb0 * bystore;
                                
                                sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) - 0.5);
                                c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                                c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                                
                                bystore = pmlArr_[kk].By_end_->point(ii,jj);
                                pmlArr_[kk].By_end_->point(ii,jj) = c_byb * pmlArr_[kk].By_end_->point(ii,jj) + c_bye * (Ez_->point((nx_-1) - ii+1, jj)-Ez_->point((nx_-1) - ii,jj));
                                Hy_->point((nx_-1) - ii,jj) = c_hyh * Hy_->point((nx_-1) - ii,jj) + c_hyb1 * pmlArr_[kk].By_end_->point(ii,jj) - c_hyb0 * bystore;
                            }
                        }
                    }
                    double eps = 1.0;
                    double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                    double sigz = 0.0;
                    double sigxx = pmlArr_[kk].sigma(0.0);
                    double sigxy = pmlArr_[kk].sigma(0.5);
                    double sigyx = 0.0;
                    double sigyy = 0.0;
                    
                    double c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                    double c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;
                    double c_hxh = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                    double c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                    double c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                    double c_byb = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                    double c_bye = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_)) ;
                    double c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                    double c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                    double c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                    if(yPML_ !=0)
                    {
                        for(int jj = a; jj < ny_ - a; jj++)
                        {
                            double bxstore = pmlArr_[kk].Bx_->point(0,jj);
                            double bystore = pmlArr_[kk].By_->point(0,jj);
                            pmlArr_[kk].Bx_->point(0,jj) = c_bxb * pmlArr_[kk].Bx_->point(0,jj) - c_bxe * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                            pmlArr_[kk].By_->point(0,jj) = c_byb * pmlArr_[kk].By_->point(0,jj) + c_bye * (Ez_->point(0+1,jj)-Ez_->point(0,jj));
                            Hx_->point(0,jj) = c_hxh * Hx_->point(0,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(0,jj) - c_hxb0 * bxstore;
                            Hy_->point(0,jj) = c_hyh * Hy_->point(0,jj) + c_hyb1 * pmlArr_[kk].By_->point(0,jj) - c_hyb0 * bystore;
                            
                            bxstore = pmlArr_[kk].Bx_end_->point(0,jj);
                            pmlArr_[kk].Bx_end_->point(0,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(0,jj) - c_bxe * (Ez_->point((nx_-1) - 0, jj+1)-Ez_->point((nx_-1) - 0,jj));
                            Hx_->point((nx_-1) - 0,jj) = c_hxh * Hx_->point((nx_-1) - 0,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(0,jj) - c_hxb0 * bxstore;
                        }
                    }
                    else
                    {
                        for(int jj = 0; jj < ny_ - 1; jj++)
                        {
                            double bxstore = pmlArr_[kk].Bx_->point(0,jj);
                            double bystore = pmlArr_[kk].By_->point(0,jj);
                            pmlArr_[kk].Bx_->point(0,jj) = c_bxb * pmlArr_[kk].Bx_->point(0,jj) - c_bxe * (Ez_->point(0,jj+1)-Ez_->point(0,jj));
                            pmlArr_[kk].By_->point(0,jj) = c_byb * pmlArr_[kk].By_->point(0,jj) + c_bye * (Ez_->point(0+1,jj)-Ez_->point(0,jj));
                            Hx_->point(0,jj) = c_hxh * Hx_->point(0,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(0,jj) - c_hxb0 * bxstore;
                            Hy_->point(0,jj) = c_hyh * Hy_->point(0,jj) + c_hyb1 * pmlArr_[kk].By_->point(0,jj) - c_hyb0 * bystore;
                            
                            bxstore = pmlArr_[kk].Bx_end_->point(0,jj);
                            pmlArr_[kk].Bx_end_->point(0,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(0,jj) - c_bxe * (Ez_->point((nx_-1) - 0, jj+1)-Ez_->point((nx_-1) - 0,jj));
                            Hx_->point((nx_-1) - 0,jj) = c_hxh * Hx_->point((nx_-1) - 0,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(0,jj) - c_hxb0 * bxstore;
                        }
                    }             
                    break;
                }
                case Y:
                {
                    if(xPML_ ==0)
                        a = 0;
                    else
                        a = xPML_;
                    double eps = 1.0;
                    double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                    double sigz = 0.0;
                    double sigxx = 0.0;
                    double sigxy = 0.0;
                    double sigyx = pmlArr_[kk].sigma(0.5);
                    double sigyy = pmlArr_[kk].sigma(0.0);
                    
                    double c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                    double c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;
                    double c_hxh = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                    double c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                    double c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                    double c_byb = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                    double c_bye = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_)) ;
                    double c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                    double c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                    double c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                     
                    if(xPML_ !=0)
                    {
                        for(int jj = a; jj < nx_ - a; jj++)
                        {
                            double bxstore = pmlArr_[kk].Bx_->point(jj,0);
                            double bystore = pmlArr_[kk].By_->point(jj,0);
                            pmlArr_[kk].Bx_->point(jj,0) = c_bxb * pmlArr_[kk].Bx_->point(jj,0) - c_bxe * (Ez_->point(jj,0+1)-Ez_->point(jj,0));
                            pmlArr_[kk].By_->point(jj,0) = c_byb * pmlArr_[kk].By_->point(jj,0) + c_bye * (Ez_->point(jj+1,0)-Ez_->point(jj,0));
                            //cout <<"H1"<<endl;
                            Hx_->point(jj,0) = c_hxh * Hx_->point(jj,0) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,0) - c_hxb0 * bxstore;
                            Hy_->point(jj,0) = c_hyh * Hy_->point(jj,0) + c_hyb1 * pmlArr_[kk].By_->point(jj,0) - c_hyb0 * bystore;
                            bystore = pmlArr_[kk].By_end_->point(jj,0);
                            pmlArr_[kk].By_end_->point(jj,0) = c_byb * pmlArr_[kk].By_end_->point(jj,0) + c_bye * (Ez_->point(jj+1, (ny_-1)-0)-Ez_->point(jj,(ny_-1)-0));
                            Hy_->point(jj,(ny_-1)-0) = c_hyh * Hy_->point(jj,(ny_-1)-0) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,0) - c_hyb0 * bystore;
                        }
                    }
                    else
                    {
                        for(int jj = 0; jj < nx_ - 1; jj++)
                        {
                            double bxstore = pmlArr_[kk].Bx_->point(jj,0);
                            double bystore = pmlArr_[kk].By_->point(jj,0);
                            pmlArr_[kk].Bx_->point(jj,0) = c_bxb * pmlArr_[kk].Bx_->point(jj,0) - c_bxe * (Ez_->point(jj,0+1)-Ez_->point(jj,0));
                            pmlArr_[kk].By_->point(jj,0) = c_byb * pmlArr_[kk].By_->point(jj,0) + c_bye * (Ez_->point(jj+1,0)-Ez_->point(jj,0));
                            //cout <<"H1"<<endl;
                            Hx_->point(jj,0) = c_hxh * Hx_->point(jj,0) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,0) - c_hxb0 * bxstore;
                            Hy_->point(jj,0) = c_hyh * Hy_->point(jj,0) + c_hyb1 * pmlArr_[kk].By_->point(jj,0) - c_hyb0 * bystore;
                            bystore = pmlArr_[kk].By_end_->point(jj,0);
                            pmlArr_[kk].By_end_->point(jj,0) = c_byb * pmlArr_[kk].By_end_->point(jj,0) + c_bye * (Ez_->point(jj+1, (ny_-1)-0)-Ez_->point(jj,(ny_-1)-0));
                            Hy_->point(jj,(ny_-1)-0) = c_hyh * Hy_->point(jj,(ny_-1)-0) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,0) - c_hyb0 * bystore;
                        }
                    }   
                    for(int ii = 1; ii < pmlArr_[kk].thickness(); ii ++)
                    {
                         eps = 1.0;
                         kapx = 1.0;  kapy = 1.0;  kapz = 1.0;
                         sigz = 0.0;
                         sigxx = 0.0;
                         sigxy = 0.0;
                         sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) + 0.5));
                         sigyy = pmlArr_[kk].sigma(static_cast<double>((ii)));
                        
                         c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                         c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;
                         c_hxh = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                         c_hxb0 = (2*eps*kapx - sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;
                         c_hxb1 = (2*eps*kapx + sigxx*dt_) / (2*eps*kapz + sigz*dt_) / eps;

                         c_byb = (2*eps*kapz - sigz*dt_) / (2*eps*kapz + sigz*dt_);
                         c_bye = 2 * eps * dt_ / (dy_ * (2*eps*kapz + sigz*dt_)) ;
                         c_hyh = (2*eps*kapx - sigxy*dt_) / (2*eps*kapx + sigxy*dt_);
                         c_hyb0 = (2*eps*kapy - sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                         c_hyb1 = (2*eps*kapy + sigyy*dt_) / (2*eps*kapx + sigxy*dt_) / eps;
                        
                        for(int jj =  a; jj < nx_ - a; jj ++)
                        {
                            if(jj != nx_-1)
                            {
                                double bxstore = pmlArr_[kk].Bx_->point(jj,ii);
                                double bystore = pmlArr_[kk].By_->point(jj,ii);
                                pmlArr_[kk].Bx_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_->point(jj,ii) - c_bxe * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                                pmlArr_[kk].By_->point(jj,ii) = c_byb * pmlArr_[kk].By_->point(jj,ii) + c_bye * (Ez_->point(jj+1,ii)-Ez_->point(jj,ii));
                                //cout <<"H1"<<endl;
                                Hx_->point(jj,ii) = c_hxh * Hx_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,ii) - c_hxb0 * bxstore;
                                Hy_->point(jj,ii) = c_hyh * Hy_->point(jj,ii) + c_hyb1 * pmlArr_[kk].By_->point(jj,ii) - c_hyb0 * bystore;

                                sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) - 0.5));
                                c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;

                                bxstore = pmlArr_[kk].Bx_end_->point(jj,ii);
                                bystore = pmlArr_[kk].By_end_->point(jj,ii);
                                pmlArr_[kk].Bx_end_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_end_->point(jj,ii) - c_bxe * (Ez_->point(jj, (ny_-1)-ii+1)-Ez_->point(jj,(ny_-1)-ii));
                                pmlArr_[kk].By_end_->point(jj,ii) = c_byb * pmlArr_[kk].By_end_->point(jj,ii) + c_bye * (Ez_->point(jj+1, (ny_-1)-ii)-Ez_->point(jj,(ny_-1)-ii));
                                //cout <<"H2"<<endl;
                                Hx_->point(jj,(ny_-1)-ii) = c_hxh * Hx_->point(jj,(ny_-1)-ii) + c_hxb1 * pmlArr_[kk].Bx_end_->point(jj,ii) - c_hxb0 * bxstore;
                                Hy_->point(jj,(ny_-1)-ii) = c_hyh * Hy_->point(jj,(ny_-1)-ii) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,ii) - c_hyb0 * bystore;
                            }
                            else
                            {
                                double bxstore = pmlArr_[kk].Bx_->point(jj,ii);

                                pmlArr_[kk].Bx_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_->point(jj,ii) - c_bxe * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                                //cout <<"H1"<<endl;
                                Hx_->point(jj,ii) = c_hxh * Hx_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,ii) - c_hxb0 * bxstore;

                                sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) - 0.5));
                                c_bxb = (2*eps*kapy - sigyx*dt_) / (2*eps*kapy + sigyx*dt_);
                                c_bxe = 2 * eps * dt_ / (dy_ * (2*eps*kapy + sigyx*dt_)) ;
                                
                                bxstore = pmlArr_[kk].Bx_end_->point(jj,ii);
                                pmlArr_[kk].Bx_end_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_end_->point(jj,ii) - c_bxe * (Ez_->point(jj, (ny_-1)-ii+1)-Ez_->point(jj,(ny_-1)-ii));
                                //cout <<"H2"<<endl;
                                Hx_->point(jj,(ny_-1)-ii) = c_hxh * Hx_->point(jj,(ny_-1)-ii) + c_hxb1 * pmlArr_[kk].Bx_end_->point(jj,ii) - c_hxb0 * bxstore;
                            }
                        }
                    }
                    
                    break;    
                }
                case Z:
                    throw logic_error("Third dimension is not implimented yet");
                    break;
                default:
                    throw logic_error("you hit the default this is bad!");
                    break;
            }
        }
    }
