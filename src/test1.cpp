for(int kk = 0; kk < pmlArr_.size(); kk++)
    {
        for(int ii = 0; ii < pmlArr_[kk].thickness(); ii ++)
        {

            switch (pmlArr_[kk].d())
            {
                case X:
                    for(int jj =  yPML_ ; jj < ny_ - yPML_; jj ++)
                    {
                        if(Ez_)
                        {
                            double eps = 1.0;
                            double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                            double sigz = 0.0;
                            double sigxx = pmlArr_[kk].sigma(static_cast<double>(ii));
                            double sigxy = pmlArr_[kk].sigma(static_cast<double>(ii) + 0.5);
                            double sigyx = 0.0;
                            double sigyy = 0.0;
                            //Kappas change throughout
                            //cout << "store1" <<endl;
                            //cout << jj << "\t" << ny_+1<< "\t"<< pmlArr_[kk].Bx_->y()  << "\t" << pmlArr_[kk].By_->y() <<endl;
                            
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
                            
                            double bxstore = pmlArr_[kk].Bx_->point(ii,jj);
                            double bystore = pmlArr_[kk].By_->point(ii,jj);
                            if(jj != ny_-1 && ii != 0)
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
                            else if(ii !=0)
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
                            else if (jj != ny_-1)
                            {
                                pmlArr_[kk].Bx_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_->point(ii,jj) - c_bxe * (Ez_->point(ii,jj+1)-Ez_->point(ii,jj));
                                pmlArr_[kk].By_->point(ii,jj) = c_byb * pmlArr_[kk].By_->point(ii,jj) + c_bye * (Ez_->point(ii+1,jj)-Ez_->point(ii,jj));
                                Hx_->point(ii,jj) = c_hxh * Hx_->point(ii,jj) + c_hxb1 * pmlArr_[kk].Bx_->point(ii,jj) - c_hxb0 * bxstore;
                                Hy_->point(ii,jj) = c_hyh * Hy_->point(ii,jj) + c_hyb1 * pmlArr_[kk].By_->point(ii,jj) - c_hyb0 * bystore;
                                bxstore = pmlArr_[kk].Bx_end_->point(ii,jj);
                                pmlArr_[kk].Bx_end_->point(ii,jj) = c_bxb * pmlArr_[kk].Bx_end_->point(ii,jj) - c_bxe * (Ez_->point((nx_-1) - ii, jj+1)-Ez_->point((nx_-1) - ii,jj));
                                Hx_->point((nx_-1) - ii,jj) = c_hxh * Hx_->point((nx_-1) - ii,jj) + c_hxb1 * pmlArr_[kk].Bx_end_->point(ii,jj) - c_hxb0 * bxstore;
                            }
                                
                        }
                        else
                        {
                            
                        }
                    }
                    break;
                case Y:
                {
                    cout << endl;
                    for(int jj = xPML_; jj < nx_ - xPML_; jj ++)
                    {
                        if(Ez_)
                        {
                            double eps = 1.0;
                            double kapx = 1.0; double kapy = 1.0; double kapz = 1.0;
                            double sigz = 0.0;
                            double sigxx = 0.0;
                            double sigxy = 0.0;
                            double sigyx = pmlArr_[kk].sigma(static_cast<double>((ii) + 0.5));
                            double sigyy = pmlArr_[kk].sigma(static_cast<double>((ii)));
                            //Kappas change throughout
                            //cout << "store1" <<endl;
                            //cout << jj << "\t" << ny_+1<< "\t"<< pmlArr_[kk].Bx_->y()  << "\t" << pmlArr_[kk].By_->y() <<endl;
                            
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

                            
                            if(ii!=0 && jj != nx_-1)
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
                            else if (jj != nx_-1)
                            {
                                double bxstore = pmlArr_[kk].Bx_->point(jj,ii);
                                double bystore = pmlArr_[kk].By_->point(jj,ii);
                                pmlArr_[kk].Bx_->point(jj,ii) = c_bxb * pmlArr_[kk].Bx_->point(jj,ii) - c_bxe * (Ez_->point(jj,ii+1)-Ez_->point(jj,ii));
                                pmlArr_[kk].By_->point(jj,ii) = c_byb * pmlArr_[kk].By_->point(jj,ii) + c_bye * (Ez_->point(jj+1,ii)-Ez_->point(jj,ii));
                                //cout <<"H1"<<endl;
                                Hx_->point(jj,ii) = c_hxh * Hx_->point(jj,ii) + c_hxb1 * pmlArr_[kk].Bx_->point(jj,ii) - c_hxb0 * bxstore;
                                Hy_->point(jj,ii) = c_hyh * Hy_->point(jj,ii) + c_hyb1 * pmlArr_[kk].By_->point(jj,ii) - c_hyb0 * bystore;
                                bystore = pmlArr_[kk].By_end_->point(jj,ii);
                                pmlArr_[kk].By_end_->point(jj,ii) = c_byb * pmlArr_[kk].By_end_->point(jj,ii) + c_bye * (Ez_->point(jj+1, (ny_-1)-ii)-Ez_->point(jj,(ny_-1)-ii));
                                Hy_->point(jj,(ny_-1)-ii) = c_hyh * Hy_->point(jj,(ny_-1)-ii) + c_hyb1 * pmlArr_[kk].By_end_->point(jj,ii) - c_hyb0 * bystore;
                            }
                            else if(ii != 0)
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
                        else
                        {
                           
                        }  
                    }
                    break;
                }
                case Z:
                    throw logic_error("While yes we could have a thrid dimension to run, I have yet to be implimented to do such a thing. So please accept this error as my sincerest appology.");
                    break;
                default:
                    throw logic_error("I would appricate it if you stick to the standard X,Y,Z directions. While it's fun to invent new ones, it is very hard to do calculations if I don't even understand what dimension I am in. Please try again!");
                    break;
            }
        }
    }