r"""
    Examples of using the class "Sbox"

AUTHORS:

- Oleksandr Kazymyrov (2011): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013 Oleksandr Kazymyrov <oleksandr.kazymyrov@ii.uib.no>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***************************************************************************** 

load ./Data.sage
load ./Sbox.sage

def test_temp(**kargs):
    bits=kargs.get('bits',3)

    S=Sbox(n=bits,m=bits)

    S.generate_sbox(method='APN6')

    print "Sbox\t\t\t\t= {0}".format(S.get_sbox())
    print "Characteristics of boolean functions:"
    print "Balanced\t\t\t= {0}".format(S.is_balanced())
    print "Nonlinearity\t\t\t= {0}".format(S.NL())
    print "Absolute indicator\t\t= {0}".format(S.absolute_indicator())
    print "Propagation criterion\t\t= {0}".format(S.PC())
    print "Correlation immunity\t\t= {0}".format(S.CI())
    print "Sum-of-squares indicator\t= {0}".format(S.SSI())
    print "Minimum degree\t\t\t= {0}".format(S.minimum_degree())
    print "Resilient\t\t\t= {0}".format(S.resilient())
    print "SAC\t\t\t\t= {0}".format(S.SAC())
    print ""
    print "Characteristics of substitution:"
    print "Bijection\t\t\t= {0}".format(S.is_bijection())
    print "Check interpolation polynomial\t= {0}".format(S.check_polynomial())
    print "Interpolation polynomial\t= {0}".format(S.interpolation_polynomial())
    print "Multiplicative generator\t= {0}".format(S.get_mg())
    print "Modulus\t\t\t\t= {0}".format(S.get_modulus())
    print "Maximum of diff table\t\t= {0}".format(S.MDT())
    print "Maximum of lin table\t\t= {0}".format(S.MLT())
    print "Cycles\t\t\t\t= {0}".format(S.cycles())
    ret=S.algebraic_immunity_sbox()
    print "Algebraic immunity\t\t: degree={0} equations={1}".format(ret[0],ret[1])
    print "Check system\t\t\t= {0}".format(S.check_system(degree=ret[0]))
    print "~"*40
    print ""

def test_all_functions(**kargs):
    bits=kargs.get('bits',3)

    Sboxes = [
        AES_sbox
    ]

    S=Sbox(n=bits,m=bits)

    for sbox in Sboxes:
        S.set_sbox(sbox)
        print "Sbox\t\t\t\t= {0}".format(S.get_sbox())
        print "Characteristics of boolean functions:"
        print "Balanced\t\t\t= {0}".format(S.is_balanced())
        print "Nonlinearity\t\t\t= {0}".format(S.NL())
        print "Absolute indicator\t\t= {0}".format(S.absolute_indicator())
        print "Propagation criterion\t\t= {0}".format(S.PC())
        print "Correlation immunity\t\t= {0}".format(S.CI())
        print "Sum-of-squares indicator\t= {0}".format(S.SSI())
        print "Minimum degree\t\t\t= {0}".format(S.minimum_degree())
        print "Resilient\t\t\t= {0}".format(S.resilient())
        print "SAC\t\t\t\t= {0}".format(S.SAC())
        print ""
        print "Characteristics of substitution:"
        print "Bijection\t\t\t= {0}".format(S.is_bijection())
        print "Check interpolation polynomial\t= {0}".format(S.check_polynomial())
        print "Interpolation polynomial\t= {0}".format(S.interpolation_polynomial())
        print "Multiplicative generator\t= {0}".format(S.get_mg())
        print "Modulus\t\t\t\t= {0}".format(S.get_modulus())
        print "Maximum of diff table\t\t= {0}".format(S.MDT())
        print "Maximum of lin table\t\t= {0}".format(S.MLT())
        print "Cycles\t\t\t\t= {0}".format(S.cycles())
        ret=S.algebraic_immunity_sbox()
        print "Algebraic immunity\t\t: degree={0} equations={1}".format(ret[0],ret[1])
        print "Check system\t\t\t= {0}".format(S.check_system(degree=ret[0]))
        print "~"*40
        print ""
                
    return 0

def test_equivalence(**kargs):
    bits=kargs.get('bits',3)

    S=Sbox(n=bits,m=bits)

    # cnt = 0

    # while(true):
    #     S.random_substitution()

    #     F=S.interpolation_polynomial()
    
    #     M1 = S.is_equivalent_to_permutation(F=F,full=True)
    #     M2 = S.is_equivalent_to_permutation_new(F=F,full=True)

    #     if len(M1) != len(M2):
    #         print "M1:\n{0}".format(M1)
    #         print "M2:\n{0}".format(M2)
    #         break            

    #     for m in M1:
    #         if not m in M2:
    #             print "M1:\n{0}".format(M1)
    #             print "M2:\n{0}".format(M2)
    #             break                

    #     cnt += 1
    #     print "cnt = {0}".format(cnt)
    # return

    #F="g*x^3+g^5*x^10+g^4*x^24" # for n=6
    F="x^3+g^4*x^10+g^3*x^24" # for n=6
    #F="g*x^3+g^257*x^514+g^256*x^528" # for n=12

    #M1 = S.is_equivalent_to_permutation(F=F,debug=True,full=True)
    #M2 = S.is_equivalent_to_permutation_new(F=F,pt=[[0,1,2,0,4,8,16,27,32,64,128,256,512,1024,267,0,0,0,0,0,0,0,0,0],[]])
    M2 = S.is_equivalent_to_permutation_new(F=F,debug=True)
    #M2 = S.is_equivalent_to_permutation_new(F=F)

    # for i,M in enumerate(M1):
    #     print "M{0}:\n{1}".format(i,M)
    # for i,M in enumerate(M2):
    #     print "M{0}:\n{1}".format(i,M)
    #print "M1:\n{0}".format(M1)
    print "M2:\n{0}".format(M2)
    #print "M2:\n{0}".format(len(M2))
    #print ">>> {0} <<<".format(M1 == M2)

    return

    S.generate_sbox(method="polynomial",G=F,T="CCZ",M=M)

    print "Sbox\t\t\t\t= {0}".format(S.get_sbox())
    print "Characteristics of boolean functions:"
    print "Balanced\t\t\t= {0}".format(S.is_balanced())
    print "Nonlinearity\t\t\t= {0}".format(S.NL())
    print "Absolute indicator\t\t= {0}".format(S.absolute_indicator())
    print "Propagation criterion\t\t= {0}".format(S.PC())
    print "Correlation immunity\t\t= {0}".format(S.CI())
    print "Sum-of-squares indicator\t= {0}".format(S.SSI())
    print "Minimum degree\t\t\t= {0}".format(S.minimum_degree())
    print "Resilient\t\t\t= {0}".format(S.resilient())
    print "SAC\t\t\t\t= {0}".format(S.SAC())
    print ""
    print "Characteristics of substitution:"
    print "Bijection\t\t\t= {0}".format(S.is_bijection())
    print "Check interpolation polynomial\t= {0}".format(S.check_polynomial())
    print "Interpolation polynomial\t= {0}".format(S.interpolation_polynomial())
    print "Multiplicative generator\t= {0}".format(S.get_mg())
    print "Modulus\t\t\t\t= {0}".format(S.get_modulus())
    print "Maximum of diff table\t\t= {0}".format(S.MDT())
    print "Maximum of lin table\t\t= {0}".format(S.MLT())
    print "Cycles\t\t\t\t= {0}".format(S.cycles())
    ret=S.algebraic_immunity_sbox()
    print "Algebraic immunity\t\t: degree={0} equations={1}".format(ret[0],ret[1])
    print "Check system\t\t\t= {0}".format(S.check_system(degree=ret[0]))
    print "~"*40
    print ""
                
    return
