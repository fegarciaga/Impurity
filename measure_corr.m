function [spin_cc, charge_cc]=measure_corr(k,l,G_up, G_up_tilde, G_dn, G_dn_tilde)
    spin_cc=0.25*(G_up(k,k)*G_up(l,l)+G_up(k,l)*G_up_tilde(k,l)-G_up(k,k)*G_dn(l,l)-G_dn(k,k)*G_up(l,l)+G_dn(k,k)*G_dn(l,l)+G_dn(k,l)*G_dn_tilde(k,l));
    spin_cc=spin_cc+0.5*(G_up(k,l)*G_dn_tilde(k,l)+G_dn(k,l)*G_up_tilde(k,l));
    charge_cc=G_up(k,k)*G_up(l,l)+G_up(k,l)*G_up_tilde(k,l)+G_up(k,k)*G_dn(l,l)+G_dn(k,k)*G_up_tilde(l,l)+G_dn(k,k)*G_dn(l,l)+G_dn(k,l)*G_dn_tilde(k,l);        
end