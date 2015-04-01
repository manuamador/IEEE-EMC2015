const c=299792458.
const mu0=4*pi*1e-7
const mur=1
const eps0=1./(mu0*c^2)
const epsr=1

using PyPlot

function urc(Q,f,f0)
    w=f/f0
    I=1/sqrt(1+Q^2*(w-1/w)^2)
    return I
end

const a=2.9#3.11
const b=3.7#3.62
const d=8.7#4.84
const f0=40e6
const f1=1000e6

const df=5e6
const nf=int((f1-f0)/df)+1
const f=linspace(f0,f1,nf)#array([20e6,40e6,45e6,56e6,65e6])#linspace(20e6,.3e9,281)#linspace(20e6,3e8,1121)
const V=a*b*d
const Q4=16*pi^2*V*f.^3/(1*c^3)
const Q1=3*V./(2*1*(2./sqrt(2*pi*f.*mur*mu0*1*1e5)))/200
const Q=1./(1./Q1+1./Q4)

const Npt=150
const x=rand(Npt)*(a-1)+0.5
const y=rand(Npt)*(b-1)+0.5
const z=rand(Npt)*(d-1)+0.5

n_pos=50
x_p=rand(n_pos)*(a-2)+1
y_p=rand(n_pos)*(b-2)+1
z_p=(rand(n_pos)*(d-2)+1)

using Grid
x_pi=InterpGrid(x_p,BCnil,InterpQuadratic)
y_pi=InterpGrid(y_p,BCnil,InterpQuadratic)
z_pi=InterpGrid(z_p,BCnil,InterpQuadratic)


z_b=z_pi[1:.2:n_pos]
x_b=x_pi[1:.2:n_pos]
y_b=y_pi[1:.2:n_pos]

n_steps=length(x_b)
figure(0)#,figsize=(8,5))
plot3D(x_b,y_b,z_b,"k.-")
grid()
axis("tight")
xlabel("\$x\$/m")
ylabel("\$y\$/m")
zlabel("\$z\$/m")
savefig("movement.pdf",bbox_inches="tight")
close()

E0=1
H0=1/(c*mu0*mur)

Nmodesx=zeros(nf)
Nx=zeros(nf)

Ex=zeros(Complex128,nf,n_steps,Npt)
for u=1:nf
  n_x=int(f[u]*(1+100/Q[u])*2*a/c)
  n_y=int(f[u]*(1+100/Q[u])*2*b/c)
  n_z=int(f[u]*(1+100/Q[u])*2*d/c)
  M=[]
  N=[]
  P=[]
  F=[]
  AMP=[]
  for mm=0:n_x
    for nn=0:n_y
      for pp=0:n_z
        fmnp=c/2*sqrt((mm/a)^2+(nn/b)^2+(pp/d)^2)
        #if abs(fmnp-f[u])<100*f[u]/Q[u]
        Amp=urc(Q[u],f[u],fmnp)
        if Amp>0.01
          M=vcat(M,mm)
          N=vcat(N,nn)
          P=vcat(P,pp)
          F=vcat(F,fmnp)
          AMP=vcat(AMP,Amp)
        end
      end
    end
  end
  nm=length(M)
  println("$(f[u]), $nm resonant freq.")
  f_coeffTMx_x=zeros(n_steps,nm)
  f_coeffTMx_y=zeros(n_steps,nm)
  f_coeffTMx_z=zeros(n_steps,nm)
  f_coeffTEx_y=zeros(Complex128,n_steps,nm)
  f_coeffTEx_z=zeros(Complex128,n_steps,nm)
  for pos=1:n_steps
    ETMx_x=zeros(Npt)
    ETMx_y=zeros(Npt)
    ETMx_z=zeros(Npt)
    ETEx_y=zeros(Complex128,Npt)
    ETEx_z=zeros(Complex128,Npt)
    xe=x_b[pos]
    ye=y_b[pos]
    ze=z_b[pos]
    #println(" Step $pos/$n_steps")
    for k=1:nm
      fmnp=F[k]
      Amp=AMP[k]
      kx=M[k]*pi/a
      ky=N[k]*pi/b
      kz=P[k]*pi/d
      k2mnp=kx^2+ky^2+kz^2
      ssc_ze=sum(sin(kx*xe).*sin(ky*ye).*cos(kz*ze))
      scs_ze=sum(sin(kx*xe).*cos(ky*ye).*sin(kz*ze))
      css_ze=sum(cos(kx*xe).*sin(ky*ye).*sin(kz*ze))
      ssc_xe=sum(sin(ky*ye).*sin(kz*ze).*cos(kx*xe))
      scs_xe=sum(sin(ky*ye).*cos(kz*ze).*sin(kx*xe))
      css_xe=sum(cos(ky*ye).*sin(kz*ze).*sin(kx*xe))
      ssc_ye=sum(sin(kz*ze).*sin(kx*xe).*cos(ky*ye))
      scs_ye=sum(sin(kz*ze).*cos(kx*xe).*sin(ky*ye))
      css_ye=sum(cos(kz*ze).*sin(kx*xe).*sin(ky*ye))

      if M[k]>0 && N[k]>0 && P[k]>0
        #TM modes z
        f_coeffTMx_z[pos,k] = -E0*Amp*kx*kz/(k2mnp-kz^2)*css_ze
        #f_coeffTMy_z[pos,k] = E0*Amp*ky*kz/(k2mnp-kz^2)*scs_ze
        #f_coeffTMz_z[pos,k] = E0*Amp*ssc_ze
        #TE modes z
        f_coeffTEx_z[pos,k] = -H0*Amp*1im*2*pi*fmnp*mur*mu0*ky/(k2mnp-kz^2)*css_ze
        #f_coeffTEy_z[pos,k] = H0*Amp*1im*2*pi*fmnp*mur*mu0*kx/(k2mnp-kz^2)*scs_ze

        #TM modes x
        #f_coeffTMy_x[pos,k] = -E0*Amp*ky*kx/(k2mnp-kx^2)*css_xe
        #f_coeffTMz_x[pos,k] = E0*Amp*kz*kx/(k2mnp-kx^2)*scs_xe
        f_coeffTMx_x[pos,k] = E0*Amp*ssc_xe
        #TE modes x
        #f_coeffTEy_x[pos,k] = -H0*Amp*1im*2*pi*fmnp*mur*mu0*kz/(k2mnp-kx^2)*css_xe
        #f_coeffTEz_x[pos,k] = H0*Amp*1im*2*pi*fmnp*mur*mu0*ky/(k2mnp-kx^2)*scs_xe

        ##TM modes y
        #f_coeffTMz_y[pos,k] += -H0*Amp*E0*Amp*kz*ky/(k2mnp-ky^2)*css_ye
        f_coeffTMx_y[pos,k] = H0*Amp*E0*Amp*kx*ky/(k2mnp-ky^2)*scs_ye
        #f_coeffTMy_y[pos,k] = H0*Amp*E0*Amp*ssc_ye
        ##TE modes y
        #f_coeffTEz_y[pos,k] = -H0*Amp*1im*2*pi*fmnp*mur*mu0*kx/(k2mnp-ky^2)*css_ye
        f_coeffTEx_y[pos,k] = H0*Amp*1im*2*pi*fmnp*mur*mu0*kz/(k2mnp-ky^2)*scs_ye
      end


      #if M[k]>0 && N[k]>0 && P[k]==0
        #f_coeffTMz_z[pos,k] = E0*Amp*ssc_ze #TM mode z
        #f_coeffTEz_y[pos,k] = -H0*Amp*1im*2*pi*fmnp*mur*mu0*kx/(k2mnp-ky^2)*css_ye #TE mode y
        #f_coeffTEz_x[pos,k] = H0*Amp*1im*2*pi*fmnp*mur*mu0*ky/(k2mnp-kx^2)*scs_xe #TE mode x
      #end

      if M[k]==0 && N[k]>0 && P[k]>0
        f_coeffTEx_z[pos,k] += -H0*Amp*1im*2*pi*fmnp*mur*mu0*ky/(k2mnp-kz^2)*css_ze #TE mode z
        f_coeffTMx_x[pos,k] += E0*Amp*ssc_xe #TM mode x
        f_coeffTEx_y[pos,k] += H0*Amp*1im*2*pi*fmnp*mur*mu0*kz/(k2mnp-ky^2)*scs_ye #TE mode y
      end

      #if M[k]>0 && N[k]==0 && P[k]>0
        #f_coeffTEy_z[pos,k] = H0*Amp*1im*2*pi*fmnp*mur*mu0*kx/(k2mnp-kz^2)*scs_ze #TE mode z
        #f_coeffTMy_y[pos,k] = E0*Amp*ssc_ye #TM mode y
        #f_coeffTEy_x[pos,k] = -H0*Amp*1im*2*pi*fmnp*mur*mu0*kz/(k2mnp-kx^2)*css_xe #TE mode x
      #end

      for i=1:Npt
        ssc_z=sin(kx*x[i])*sin(ky*y[i])*cos(kz*z[i])
        scs_z=sin(kx*x[i])*cos(ky*y[i])*sin(kz*z[i])
        css_z=cos(kx*x[i])*sin(ky*y[i])*sin(kz*z[i])
        ssc_x=sin(ky*y[i])*sin(kz*z[i])*cos(kx*x[i])
        scs_x=sin(ky*y[i])*cos(kz*z[i])*sin(kx*x[i])
        css_x=cos(ky*y[i])*sin(kz*z[i])*sin(kx*x[i])
        ssc_y=sin(kz*z[i])*sin(kx*x[i])*cos(ky*y[i])
        scs_y=sin(kz*z[i])*cos(kx*x[i])*sin(ky*y[i])
        css_y=cos(kz*z[i])*sin(kx*x[i])*sin(ky*y[i])

        if M[k]>0 && N[k]>0 && P[k]>0#TM modes z
          #TM modes z
          ETMx_z[i] += f_coeffTMx_z[pos,k]*css_z
          #ETMy_z[i] += f_coeffTMy_z[pos,k]*scs_z
          #ETMz_z[i] += f_coeffTMz_z[pos,k]*ssc_z

          #TE modes z
          ETEx_z[i] += f_coeffTEx_z[pos,k]*css_z
          #ETEy_z[i] += f_coeffTEy_z[pos,k]*scs_z

          #TM modes x
          #ETMy_x[i] += f_coeffTMy_x[pos,k]*css_x
          #ETMz_x[i] += f_coeffTMz_x[pos,k]*scs_x
          ETMx_x[i] += f_coeffTMx_x[pos,k]*ssc_x

          #TE modes x
          #ETEy_x[i] += f_coeffTEy_x[pos,k]*css_x
          #ETEz_x[i] += f_coeffTEz_x[pos,k]*scs_x

          ##TM modes y
          #ETMz_y[i] += f_coeffTMz_y[pos,k]*css_y
          ETMx_y[i] += f_coeffTMx_y[pos,k]*scs_y
          #ETMy_y[i] += f_coeffTMy_y[pos,k]*ssc_y

          ##TE modes y
          #ETEz_y[i] += f_coeffTEz_y[pos,k]*css_y
          ETEx_y[i] += f_coeffTEx_y[pos,k]*scs_y
        end

        #if M[k]>0 && N[k]>0 && P[k]==0
          #ETMz_z[i] += f_coeffTMz_z[pos,k]*ssc_z #TM mode z
          #ETEz_y[i]+= f_coeffTEz_y[pos,k]*css_y #TE mode y
          #ETEz_x[i]+= f_coeffTEz_x[pos,k]*scs_x #TE mode x
        #end

        if M[k]==0 && N[k]>0 && P[k]>0
          ETEx_z[i] += f_coeffTEx_z[pos,k]*css_z #TE mode z
          ETMx_x[i] += f_coeffTMx_x[pos,k]*ssc_x #TM mode x
          ETEx_y[i] += f_coeffTEx_y[pos,k]*scs_y #TE mode y
        end

        #if M[k]>0 && N[k]==0 && P[k]>0
          #ETEy_z[i] += f_coeffTEy_z[pos,k]*scs_z #TE mode z
          #ETMy_y[i] += f_coeffTMy_y[pos,k]*ssc_y #TM mode y
          #ETEy_x[i] += f_coeffTEy_x[pos,k]*css_x #TE mode x
        #end
      end
    end
    Ex[u,pos,:]=ETMx_x+ETMx_y+ETMx_z+ETEx_y+ETEx_z
    #Ey[u,pos,:]=ETMy_x+ETMy_y+ETMy_z+ETEy_x+ETEy_z
    #Ez[u,pos,:]=ETMz_x+ETMz_y+ETMz_z+ETEz_y+ETEz_x
  end


  f_coeffTMx=f_coeffTMx_x+f_coeffTMx_y+f_coeffTMx_z
  f_coeffTEx=f_coeffTEx_y+f_coeffTEx_z
  modesTE=zeros(n_z)
  modesTM=zeros(n_z)
  maxmax=maximum([maximum(abs(f_coeffTMx)),maximum(abs(f_coeffTEx))])
  for k=1:nm
    if maximum(abs(f_coeffTEx[:,k]))/maxmax>0.3
      modesTE[P[k]] += 1
    end
    if maximum(abs(f_coeffTMx[:,k]))/maxmax>0.3
      modesTM[P[k]] += 1
    end
  end


  ax=cor(squeeze(Ex[u,:,:],1)')
  Nx[u]=n_steps^2/sum(abs(ax).^2)
  Nmodesx[u]=sum((modesTE+modesTM).>0)
  println("Nx=$(Nx[u]), Nmodesx=$(Nmodesx[u])")


  figure(num=1,figsize=(14,10))
  ftext=round(f[u]*1000/1e6)/1000
  suptitle("\$ f=$ftext\$ MHz",fontsize=15)
  subplot(211)
  Nxtxt=round(Nx[u]*100)/100
  title("\$ \\vert E_x\\vert\$, \$ N_x=$Nxtxt \$ pos.")
  #imshow(abs(ax))
  #xlabel("\$z\$")
  #ylabel("\$z\$/cm")
  #colorbar(shrink=0.7)
  #clim(0,1)
  plot(abs(squeeze(Ex[u,:,:],1)),"k",lw=0.3)

  xlabel("steps")
  ylabel("V/m")
  grid()

  subplot(212)
  if Nmodesx[u]>1
    title("Main modes, \$N_m=$(Nmodesx[u])\$ uncorrelated modes")
  else
    title("Main modes, \$N_m=$(Nmodesx[u])\$ uncorrelated mode")
  end

  for i=1:nm
    if maximum(abs(f_coeffTEx[:,i]))/maxmax.>0.3
      Mtxt=M[i]
      Ntxt=N[i]
      Ptxt=P[i]
      label="$Mtxt \\ $Ntxt \\ $Ptxt"
      plot(imag(f_coeffTEx[:,i])/maxmax,label="\$ TE^{yz}_{$label}\$")
    else
      plot(imag(f_coeffTEx[:,i])/maxmax,"k",lw=0.3)
    end

    if maximum(abs(f_coeffTMx[:,i]))/maxmax.>0.3
      Mtxt=M[i]
      Ntxt=N[i]
      Ptxt=P[i]
      label="$Mtxt \\ $Ntxt \\ $Ptxt"
      plot(f_coeffTMx[:,i]/maxmax,"--",label="\$ TM^{xyz}_{$label}\$")
    else
      plot(f_coeffTMx[:,i]/maxmax,"--k",lw=0.3)
    end
  end

  ylim(-1.1,1.1)
  grid()
  xlabel("steps")
  ylabel("A.U.")
  legend(fontsize="small")
  savefig("N_1sourcerand$u.png")
  clf()
  M=[]
  N=[]
  P=[]
end

plot(f,Nx)
plot(f,Nmodesx)
grid()

using NPZ
npzwrite("Exbrownian1s5.npz", ["Ex" => Ex, "f" => f, "Q" => Q, "x_b" => x_b,"y_b" => y_b,"z_b" => z_b,"x" => x,"y" => y,"z" => z,"Nx"=>Nx,"Nmodesx"=>Nmodesx])
