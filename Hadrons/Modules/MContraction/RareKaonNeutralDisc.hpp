/*
 * RareKaonNeutralDisc.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */

#ifndef Hadrons_MContraction_RareKaonNeutralDisc_hpp_
#define Hadrons_MContraction_RareKaonNeutralDisc_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Weak Hamiltonian + current contractions, disconnected topology for neutral 
 * mesons.
 * 
 * These contractions are generated by operators Q_1,...,10 of the dS=1 Weak
 * Hamiltonian in the physical basis and an additional current J (see e.g. 
 * Fig 11 of arXiv:1507.03094).
 * 
 * Schematic:
 *                        
 *           q2          q4                q3
 *       /---<--¬     /---<--¬          /---<--¬
 *      /        \   /        \        /        \
 *     /          \ /          \      /          \
 *  i *            * G          |  J *            * f
 *     \           * G         /      \          /
 *      \         / \         /        \        /
 *       \--->---/   \-------/          \------/
 *          q1 
 *                      two traces
 * ----------------------------------------------------
 *
 *           q2          q4                q3
 *       /---<--¬      /---<--¬          /---<--¬
 *      /        \    /        \        /        \
 *     /          \  /          \      /          \
 *  i *         G *  * G         |  J *            * f
 *     \           /\           /      \          /
 *      \         /  \         /        \        /
 *       \--->---/    \-------/          \------/
 *          q1 
 *                      three traces
 * options
 * - q1: input propagator 1 (string)
 * - q2: input propagator 2 (string)
 * - q3: input propagator 3 (string), assumed to be sequential propagator 
 * - q4: input propagator 4 (string), assumed to be a loop
 * 
 * two traces:   trace(q1*adj(q2)*g5*G*loop*G)*trace(q3*g5)
 * three traces: trace(q1*adj(q2)*g5*G)*trace(loop*G)*trace(q3*g5)
 */

/******************************************************************************
 *                         RareKaonNeutralDisc                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

// Neutral 4pt disconnected subdiagram contractions.
#define MAKE_DISC_MESON(Q_1, Q_2, gamma) (Q_1*adj(Q_2)*g5*gamma)
#define MAKE_DISC_LOOP(Q_LOOP, gamma) (Q_LOOP*gamma)
#define MAKE_DISC_CURR(Q_c, gamma) (trace(Q_c*gamma))
//// Sum and store correlator.
#define MAKE_DIAG(exp, buf, res, n)\
sliceSum(exp, buf, Tp);\
res.name = (n);\
res.corr.resize(buf.size());\
for (unsigned int t = 0; t < buf.size(); ++t)\
{\
    res.corr[t] = TensorRemove(buf[t]);\
}

//// Contraction of mu index: use 'mu' variable in exp.
#define SUM_MU(buf,exp)\
   buf = Zero();                \
for (unsigned int mu = 0; mu < ndim; ++mu)\
{\
    buf += exp;\
}


class RareKaonNeutralDiscPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RareKaonNeutralDiscPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, q3,
                                    std::string, q4,
                                    std::string, output);
};

template <typename FImpl>
class TRareKaonNeutralDisc: public Module<RareKaonNeutralDiscPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, op,
                                        unsigned int,   trace);
    };
    typedef Correlator<Metadata> Result;
public:
    /* constructor */ 
    TRareKaonNeutralDisc(const std::string name);
    /* destructor */ 
    virtual ~TRareKaonNeutralDisc(void) {};
    /* dependency relation */ 
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    /* setup */ 
    virtual void setup(void);
    /* execution */ 
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RareKaonNeutralDisc, TRareKaonNeutralDisc<FIMPL>, MContraction);

/*******************************************************************************
 *                  TRareKaonNeutralDisc implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TRareKaonNeutralDisc<FImpl>::TRareKaonNeutralDisc(const std::string name)
: Module<RareKaonNeutralDiscPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TRareKaonNeutralDisc<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q1, par().q2,
                                   par().q3, par().q4};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TRareKaonNeutralDisc<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRareKaonNeutralDisc<FImpl>::setup(void)
{
    envTmpLat(ComplexField, "corr");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRareKaonNeutralDisc<FImpl>::execute(void)
{
    LOG(Message) << "Computing Weak Hamiltonian neutral disconnected contractions '" 
                 << getName() << "' using quarks '" << par().q1 << "', '" 
                 << par().q2 << ", '" << par().q3 << "' and '" << par().q4 
                 << "'." << std::endl;

    std::vector<Result>   result;
    Result                r;
    auto                  &q1 = envGet(PropagatorField, par().q1);
    auto                  &q2 = envGet(PropagatorField, par().q2);
    auto                  &q3 = envGet(PropagatorField, par().q3);
    auto                  &q4 = envGet(PropagatorField, par().q4);
    Gamma                 g5(Gamma::Algebra::Gamma5);

    envGetTmp(ComplexField, corr);
    for (auto &G: Gamma::gall)
    {
        SlicedComplex buf;

        r.info.op = G.g;
        // two traces
        corr = trace(q1*adj(q2)*g5*G*q4*G)*trace(q3*g5);
        sliceSum(corr, buf, Tp);
        r.corr.clear();
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            r.corr.push_back(TensorRemove(buf[t]));
        }
        r.info.trace = 2;
        result.push_back(r);
        // three traces
        corr = trace(q1*adj(q2)*g5*G)*trace(q4*G)*trace(q3*g5);
        sliceSum(corr, buf, Tp);
        r.corr.clear();
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            r.corr.push_back(TensorRemove(buf[t]));
        }
        r.info.trace = 3;
        result.push_back(r);
    }
    // IO
    saveResult(par().output, "RK_disc0", result);
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_RareKaonNeutralDisc_hpp_
