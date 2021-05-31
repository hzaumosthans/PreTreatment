import org.openscience.cdk.*;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.xmlcml.euclid.Int;

import javax.naming.InvalidNameException;
import java.time.temporal.ChronoField;
import java.util.List;
import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class Mymol {
    public IAtomContainer mol;
    private String smiles;
    private String outfileroad;
    private String sdffile;
    private FileWriter fw = null;
    private PrintWriter pw;
    private FileWriter fw2 = null;
    private PrintWriter pw2;
    private int AtomCount;
    private boolean IsUSRIn;
    private String[] sms;
    //private IAtom[] Atoms = new Atom[AtomCount];

    public Mymol(IAtomContainer molin, String of,String sf,boolean IsUSRIn){
        mol = molin;
        AtomCount = mol.getAtomCount();
        smiles = GetSmiles(mol);
        outfileroad = of;
        sdffile = sf;
        this.IsUSRIn = IsUSRIn;
        try {
            //如果文件存在，则追加内容；如果文件不存在，则创建文件
            File f=new File(outfileroad);
            fw = new FileWriter(f, true);
        } catch (IOException e) {
            e.printStackTrace();
        }
        pw = new PrintWriter(fw);
        try {
            //如果文件存在，则追加内容；如果文件不存在，则创建文件
            File f2=new File(sdffile);
            fw2 = new FileWriter(f2, true);
        } catch (IOException e) {
            e.printStackTrace();
        }
        pw2 = new PrintWriter(fw2);
    }

    String returnsmiles(){
        return smiles;
    }

    void InPutSdf(String smiles){
        if(!smiles.equals("")){
            try {
                SmilesParser sp2  = new SmilesParser(SilentChemObjectBuilder.getInstance());
                IAtomContainer m2   = sp2.parseSmiles(smiles);
                BufferedWriter bw2 = new BufferedWriter(fw2);
                SDFWriter sdfw2 = new SDFWriter(bw2);
                sdfw2.write(m2);
                bw2.flush();
            } catch (InvalidSmilesException e) {
                System.err.println(e.getMessage());
            } catch (CDKException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }


    }

    String GetSmiles(IAtomContainer mol){
        SmilesGenerator generator = SmilesGenerator.generic();
        String smiles = generator.createSMILES(mol);
        return smiles;
    }

    public void BondCountTest(){

        Map<String, Integer> bcmap = new HashMap<String, Integer>();
        bcmap.put("C",4);
        bcmap.put("H",1);
        bcmap.put("O",2);
        bcmap.put("S",2);
        bcmap.put("N",3);
        bcmap.put("Cl",1);
        bcmap.put("F",1);
        bcmap.put("Br",1);
        bcmap.put("I",1);
        bcmap.put("P",5);
        bcmap.put("Si",4);
        bcmap.put("Na",1);
        bcmap.put("Pd",7);
        bcmap.put("B",3);
        bcmap.put("Sc",8);
        bcmap.put("V",8);
        bcmap.put("Ge",8);
        int AtomId = 0;
        for (IAtom atom : mol.atoms()) {
            int BondCount = 0;
            for(IBond bond : mol.getConnectedBondsList(atom)){
                switch (bond.getOrder()) {
                    case DOUBLE -> BondCount += 2;
                    case SINGLE -> BondCount += 1;
                    case TRIPLE -> BondCount += 3;
                    case QUADRUPLE -> BondCount += 4;
                }
            }
            //System.out.println(atom.getSymbol());
            if(BondCount>bcmap.get(atom.getSymbol())){
                //System.out.println("Atom "+AtomId+" : "+"The "+atom.getSymbol()+" has "+BondCount+" bonds.");
                pw.println("Availability of molecule: ID is "+AtomId+" : "+"The "+atom.getSymbol()+" has "+BondCount+" bonds.");
                pw.flush();
                try {
                    fw.flush();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            AtomId += 1;
        }
        //System.out.print("Number of hydrogens:"+hydrogenCount);
        //System.out.print(mol.getConnectedAtomsCount(mol.getAtom(0)));


    }

    public boolean SpecialSelect(){
        int AtomId = 0;
        boolean ret = false;
        int Bcount = 0;
        for (IAtom atom : mol.atoms()) {
            if (atom.getSymbol().equals("Ge")) {
                pw.println("Availability of molecule: " + AtomId + " is " + atom.getSymbol() + ".");
                pw.flush();
                ret = true;
            }
            if (atom.getSymbol().equals("B")) {
                Bcount++;
            }

        }

        if(Bcount>=6){
            pw.println("Availability of molecule: "+Bcount+" B in the molecule.");
            pw.flush();
            ret = true;
        }
        return ret;

    }

    public void TautomerTest(){
        int AtomId = 0;
        for (IAtom atom : mol.atoms()) {
            int BondCount = 0;
            for(IBond bond : mol.getConnectedBondsList(atom)){
                switch (bond.getOrder()) {
                    case DOUBLE -> BondCount = 2;
                    case SINGLE -> BondCount = 1;
                    case TRIPLE -> BondCount = 3;
                    case QUADRUPLE -> BondCount = 4;
                }
                if(BondCount==2){

                    IAtom NextA = bond.getBegin();
                    IAtom NextB = bond.getEnd();
                    if(NextA.getSymbol().equals("O")){
                        if(NextB.getSymbol().equals("C")){
                            int hc = NextB.getImplicitHydrogenCount();
                            if(hc>0){
                                //System.out.println("A");
                                pw.println("Tautomer enumeration: ID : "+AtomId+ " Symbol: " +atom.getSymbol());
                                pw.flush();
                                try {
                                    fw.flush();
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }
                                continue;
                            }
                            for(IAtom ia : mol.getConnectedAtomsList(NextB)){
                                int bc = 0;
                                for(IBond bonds : mol.getConnectedBondsList(NextB)) {
                                    switch (bonds.getOrder()) {
                                        case DOUBLE -> bc += 2;
                                        case SINGLE -> bc += 1;
                                        case TRIPLE -> bc += 3;
                                        case QUADRUPLE -> bc += 4;
                                    }
                                }
                                if(ia.getSymbol().equals("H")){
                                    pw.println("Atom "+AtomId+ " " + atom.getSymbol() + " : "+"It is a tautomer.");
                                    pw.flush();
                                    try {
                                        fw.flush();
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }
                                    break;
                                }
                                else if(bc<4){
                                    pw.println("Atom "+AtomId+ " " +atom.getSymbol() + " : "+"It is a tautomer.");
                                    pw.flush();
                                    try {
                                        fw.flush();
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }
                                    break;
                                }
                            }

                        }
                        else{
                            continue;
                        }
                    }

                    if(NextB.getSymbol().equals("O")){
                        if(NextA.getSymbol().equals("C")){
                            int hc = NextB.getImplicitHydrogenCount();
                            if(hc>0){
                                pw.println("Atom "+AtomId+ " " +atom.getSymbol() + " : "+"It is a tautomer.");
                                pw.flush();
                                try {
                                    fw.flush();
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }
                                continue;
                            }
                            for(IAtom ia : mol.getConnectedAtomsList(NextA)){
                                int bc = 0;
                                for(IBond bonds : mol.getConnectedBondsList(NextA)) {
                                    switch (bonds.getOrder()) {
                                        case DOUBLE -> bc += 2;
                                        case SINGLE -> bc += 1;
                                        case TRIPLE -> bc += 3;
                                        case QUADRUPLE -> bc += 4;
                                    }
                                }
                                if(ia.getSymbol().equals("H")){
                                    pw.println("Atom "+AtomId+ " " +atom.getSymbol() + " : "+"It is a tautomer.");
                                    pw.flush();
                                    try {
                                        fw.flush();
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }
                                    break;
                                }
                                else if(bc<4){
                                    pw.println("Atom "+AtomId+ " " +atom.getSymbol() + " : "+"It is a tautomer.");
                                    pw.flush();
                                    try {
                                        fw.flush();
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }
                                    break;
                                }
                            }
                        }
                        else{
                            continue;
                        }
                    }


                }
            }
            //System.out.println(atom.getSymbol());
            AtomId += 1;
        }
    }

    public boolean HandTest(int iA,IAtom IA){
        mol.getConnectedBondsCount(iA);
        if(mol.getConnectedBondsCount(IA)<2){
            return false;
        }
        if(mol.getConnectedBondsCount(IA)==3){
            List<IAtom> conlist0 = mol.getConnectedAtomsList(IA);
            if(conlist0.get(0).getSymbol().equals(conlist0.get(1).getSymbol())&&!conlist0.get(0).getSymbol().equals('C')){
                return false;
            }
            if(conlist0.get(0).getSymbol().equals(conlist0.get(2).getSymbol())&&!conlist0.get(0).getSymbol().equals('C')){
                return false;
            }
            if(conlist0.get(1).getSymbol().equals(conlist0.get(2).getSymbol())&&!conlist0.get(0).getSymbol().equals('C')){
                return false;
            }

            int restatomcount = AtomCount;
            int conatomnum[] = new int[AtomCount];
            for(int i:conatomnum){
                i = -1;
            }
            conatomnum[iA] = 0;
            for(int i = 0;i<3;i++){
                List<IAtom> conlist = mol.getConnectedAtomsList(IA);

                conatomnum[mol.getAtomNumber(conlist.get(i))] = i+1;
            }
            restatomcount -= 4;
            while(restatomcount>0){
                for(int i = 0;i<AtomCount;i++){
                    if(conatomnum[i]!=-1){

                        List<IAtom> conlist = mol.getConnectedAtomsList(mol.getAtom(i));
                        int thisconnectnum = mol.getConnectedBondsCount(mol.getAtom(i));
                        for(int j = 0; j<thisconnectnum;j++){
                            if(conatomnum[conlist.get(j).getIndex()]==-1){
                                conatomnum[conlist.get(j).getIndex()] = conatomnum[i];
                                restatomcount -= 1;
                            }
                        }
                    }
                }
            }
            double sum1 = 0.0;
            double sum2 = 0.0;
            double sum3 = 0.0;
            for(int i = 0;i<AtomCount;i++){
                if(conatomnum[i]==1){sum1 += mol.getAtom(i).getMassNumber();}
                if(conatomnum[i]==2){sum2 += mol.getAtom(i).getMassNumber();}
                if(conatomnum[i]==3){sum3 += mol.getAtom(i).getMassNumber();}

            }

            if(sum1!=sum2&&sum2!=sum3&&sum1!=sum3){
                return true;
            }
            else return false;


        }
        if(mol.getConnectedBondsCount(IA)==4){
            List<IAtom> conlist0 = mol.getConnectedAtomsList(IA);

            int restatomcount = AtomCount;
            int conatomnum[] = new int[AtomCount];
            for(int i:conatomnum){
                i = -1;
            }
            conatomnum[iA] = 0;
            for(int i = 0;i<4;i++){
                List<IAtom> conlist = mol.getConnectedAtomsList(IA);
                conatomnum[mol.getAtomNumber(conlist.get(i))] = i+1;
            }
            restatomcount -= 5;
            while(restatomcount>0){
                for(int i = 0;i<AtomCount;i++){
                    if(conatomnum[i]!=-1){

                        List<IAtom> conlist = mol.getConnectedAtomsList(mol.getAtom(i));
                        int thisconnectnum = mol.getConnectedBondsCount(mol.getAtom(i));
                        for(int j = 0; j<thisconnectnum;j++){
                            if(conatomnum[mol.getAtomNumber(conlist.get(i))]==-1){
                                conatomnum[mol.getAtomNumber(conlist.get(i))] = conatomnum[i];
                                restatomcount -= 1;
                            }
                        }
                    }
                }
            }
            double sum1 = 0.0;
            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            for(int i = 0;i<AtomCount;i++){
                if(conatomnum[i]==1){sum1 += mol.getAtom(i).getMassNumber();}
                if(conatomnum[i]==2){sum2 += mol.getAtom(i).getMassNumber();}
                if(conatomnum[i]==3){sum3 += mol.getAtom(i).getMassNumber();}
                if(conatomnum[i]==4){sum4 += mol.getAtom(i).getMassNumber();}
            }
            if(sum1!=sum2&&sum2!=sum3&&sum1!=sum3||sum1!=sum3&&sum4!=sum3&&sum1!=sum4||sum3!=sum2&&sum2!=sum4&&sum4!=sum3
                    ||sum1!=sum2&&sum2!=sum4&&sum1!=sum4){
                return true;
            }
            else return false;


        }
        System.out.println("Some unexcepted problem.");
        return false;


    }

    public String DelHydro(String smiles){
        String tmp = "";
        String p = "";
        for(int i = 0;i<smiles.length();i++){
            Character tmps = smiles.charAt(i);
            if(tmps.equals('H')||tmps.equals('h')){
                continue;
            }
            tmp += tmps;
        }
        for(int i = 0;i<tmp.length();i++){
            Character tmps1 = tmp.charAt(i);
            Character tmps2;
            if(i<tmp.length()-1){
                tmps2 = tmp.charAt(i+1);
            }
            else{
                tmps2 = '*';
            }
            if(tmps1.equals('(')&&tmps2.equals(')')){
                i++;
                continue;
            }
            p += tmps1;

        }

        return p;
    }

    public boolean IsASRING(String smiles){
        boolean flag = true;
        int m = 0;
        int n = 0;
        Character cx;
        Character cy = smiles.charAt(0);
        for(m=0;m<smiles.length();m++){
            cx = smiles.charAt(m);
            if(Character.isDigit(cx)){
                for(n = 0;n<smiles.length();n++){
                    if(m==n){
                        continue;
                    }
                    cy = smiles.charAt(n);
                    if(cx.equals(cy)){
                        break;
                    }
                }
                //System.out.println("m: "+m+"  n:"+n+"  l:"+smiles.length());
                if(n==smiles.length()&&!cx.equals(cy)){
                    flag = false;
                }
            }
            if(flag==false){
                break;
            }
        }
        return flag;
    }

    public String RingStandard(String smiles){
        String tmpsmiles = "";
        int count = 1;
        int[] marks = new int[200];
        for(int i = 0;i<200;i++){
            marks[i] = 0;
        }
        int markcount = 0;
        int marklen = 0;
        for(int i = 0;i<smiles.length();i++) {
            if (Character.isDigit(smiles.charAt(i))) {
                int m = 0;
                for (int j = 0; j < marklen; j++) {
                    if (marks[m] == i) {
                        break;
                    }
                    m++;
                }
                if (m != marklen) {
                    continue;
                } else {
                    Character tmpmark = smiles.charAt(i);
                    int j = i + 1;
                    while (j<smiles.length()&&(!Character.isDigit(smiles.charAt(j)) || !tmpmark.equals(smiles.charAt(j))))  {
                        j++;
                    }
                    marks[2 * markcount] = i;
                    marks[2 * markcount + 1] = j;
                    markcount += 1;
                    marklen += 2;
                }
            }
        }
        for(int i = 0;i<smiles.length();i++){
            if (!Character.isDigit(smiles.charAt(i))) {
                tmpsmiles += smiles.charAt(i);
            }
            else{
                int marknum = 0;
                for(int j = 0; j < marklen; j++){
                    if(i==marks[j]){
                        if(j%2==0){
                            tmpsmiles += count;
                            count += 1;
                        }
                        else{
                            tmpsmiles += tmpsmiles.charAt(marks[j-1]);
                        }
                    }
                }

            }
        }



        return tmpsmiles;
    }

    public int[] SubBranchFind(String smiles,int startpos){
        //返回的是（的位置
        if(startpos>=smiles.length()-1){
            int[] last = {0};
            return last;
        }
        Character tmp = smiles.charAt(startpos);
        if(!tmp.equals('(')){
            int[] last = {0};
            return last;
        }

        int pos[] = new int[100];
        for(int m = 0;m<100;m++){
            pos[m] = 0;
        }
        int leftcount = 1;
        int rightcount = 0;
        int i = startpos + 2;
        int count = 0;
        int tmppos = 0;
        while(leftcount!=rightcount){
            Character tmpc = smiles.charAt(i);

            if(tmpc.equals('(')){
                leftcount++;

            }
            if(tmpc.equals(')')) {
                rightcount++;
                if(leftcount==rightcount){
                    pos[count] = startpos+1;
                    count++;
                    tmppos = i+1;
                    i++;

                    break;
                }
            }

            i++;
        }
        //System.out.println(smiles);
        Character tmpc2 = smiles.charAt(i);

        //System.out.println("huohuo:"+i+tmpc2);
        if(tmpc2.equals('(')){
            i++;

            leftcount = 1;
            rightcount = 0;
            while(leftcount!=rightcount){
                Character tmpc = smiles.charAt(i);
                if(tmpc.equals('(')){
                    leftcount++;

                }
                if(tmpc.equals(')')) {
                    rightcount++;
                    if(leftcount==rightcount){
                        pos[count] = tmppos;
                        count++;
                        i++;
                        break;
                    }
                }
                i++;
            }
            Character ctmp = smiles.charAt(smiles.length()-1);
            if(!ctmp.equals(')')){
                pos[count] = i;

            }
            else {count--;}


        }
        if(tmpc2.equals(')')){
            int[] last = {-1};
            return last;
        }
        else{
            /*
            if(i!=smiles.length()-1){
                pos[count] = i;
            }
            else{count--;}

             */
            pos[count] = i;
        }

        if(count>0){
            int[] last = new int[count+1];
            for(int m = 0;m<last.length;m++){
                last[m] = pos[m];
                //System.out.println(smiles+" "+m+" : "+last[m]+" start: "+startpos);
            }
            return last;
        }
        else{
            int[] last = {0};
            return last;
        }
    }

    public String CatchBranch(String smiles,int startpos){
        //startpos从（开始
        int leftcount = 1;
        int rightcount = 0;
        int i = startpos+1;
        String tmpsmiles = "";
        if(startpos==smiles.length()-1){
            return "" + smiles.charAt(i-1);

        }
        Character tmps2 = smiles.charAt(startpos);
        Character tmps;
        if(!tmps2.equals('(')&&!tmps2.equals(')')){
            int balance = 1;
            int count = startpos;
            while(balance>0&&count<smiles.length()){
                Character cx = smiles.charAt(count);
                if(cx.equals('(')){balance++;}
                if(cx.equals(')')){balance--;}
                if(balance>0){tmpsmiles += cx;}
                count++;
            }

            return tmpsmiles;

        }
        if(tmps2.equals(')')){
            return "";
        }
        tmps = smiles.charAt(i);
        while(leftcount!=rightcount&&smiles.length()>i){
            tmps = smiles.charAt(i);
            if(tmps.equals('(')){
                leftcount++;
            }
            if(tmps.equals(')')){
                rightcount++;
            }
            if(leftcount!=rightcount) {
                tmpsmiles += tmps;
            }
            i++;
        }
        return  tmpsmiles;
    }

    public String ReservePreBranch(String smiles,int preendpos){
        String tmpsmiles = "";
        Character tmpc1;
        Character tmpc2;
        Character tmpc3;
        int[] bpos;
        int i = 0;
        int flagpos = 0;
        tmpc1 = smiles.charAt(preendpos);
        if(preendpos>smiles.length()-1){return "";}
        if(tmpc1.equals('(')||tmpc1.equals('[')){
            if(preendpos>0){
                tmpc2 = smiles.charAt(preendpos-1);
                if(Character.isDigit(tmpc2)){
                    flagpos = 2;
                }
                else{
                    flagpos = 1;
                }
            }
            else {
                flagpos = 0;
            }
        }
        else{flagpos = 0;}
        //System.out.println("flagpos");
        while(i<preendpos-flagpos){
            //System.out.println(i);
            tmpc1 = smiles.charAt(i);
            if(i+2==preendpos-flagpos){
                tmpc2 = smiles.charAt(i+1);
                if(tmpc2>'a'&&tmpc2<'z'||Character.isDigit(tmpc2)){
                    tmpsmiles = "" + tmpc2 + tmpc1 + tmpsmiles;
                    return tmpsmiles;
                }
                else{
                    tmpsmiles = "" + tmpc1 + tmpc2 + tmpsmiles;
                    return tmpsmiles;
                }
            }
            if(i+1==preendpos-flagpos){
                tmpsmiles = "" + tmpc1 + tmpsmiles;
                return tmpsmiles;
            }
            if(i+1<preendpos-flagpos){
                tmpc2 = smiles.charAt(i+1);
                tmpc3 = smiles.charAt(i+2);
                if(tmpc3.equals('(')){
                    bpos = SubBranchFind(smiles,i+2);
                    if(bpos.length==1){
                        tmpsmiles = "" + tmpc1 + tmpsmiles;
                        i++;
                        continue;
                    }
                    if(bpos.length==2){
                        if(tmpc2>'a'&&tmpc2<'z'||Character.isDigit(tmpc2)){
                            tmpsmiles = "" + tmpc1 + tmpc2 + "(" + CatchBranch(smiles,bpos[0]) + ")" + tmpsmiles;
                            i = bpos[1];
                            continue;
                        }
                        else{
                            tmpsmiles = "" + tmpc2 + "(" + CatchBranch(smiles,bpos[0]) + ")" + tmpc1 + tmpsmiles;
                            i = bpos[1];
                            continue;
                        }

                    }
                    if(bpos.length==3){
                        if(tmpc2>'a'&&tmpc2<'z'||Character.isDigit(tmpc2)){
                            tmpsmiles = "" + tmpc1 + tmpc2 + "(" + CatchBranch(smiles,bpos[0]) + ")" + "(" + CatchBranch(smiles,bpos[1]) + ")" + tmpsmiles;
                            i = bpos[2];
                            continue;
                        }
                        else{
                            tmpsmiles = "" + tmpc2 + "(" + CatchBranch(smiles,bpos[0]) + ")" + "(" + CatchBranch(smiles,bpos[1]) + ")" + tmpc1 + tmpsmiles;
                            i = bpos[2];
                            continue;
                        }
                    }
                }
                if(tmpc1.equals('[')){

                    int j = i;
                    while(!tmpc1.equals(']')){
                        i++;
                        tmpc1 = smiles.charAt(i);
                    }
                    i++;
                    tmpc1 = smiles.charAt(i);
                    if(tmpc1.equals('[')){
                        int k = i;
                        while(!tmpc1.equals(']')){
                            i++;
                            tmpc1 = smiles.charAt(i);
                        }
                        if(j==0){
                            tmpsmiles = smiles.substring(k,i+1) + smiles.substring(j,k);
                            i++;
                            continue;
                        }
                        if(j==1){
                            tmpsmiles = smiles.substring(k,i+1) + smiles.substring(j,k) + smiles.substring(0,j);
                            i++;
                            continue;
                        }
                        if(j>1){
                            Character cx = tmpsmiles.charAt(j-1);
                            if(cx>'a'&&cx<'z'||Character.isDigit(cx)){
                                tmpsmiles = tmpsmiles.substring(0,2)+smiles.substring(k,i+1) + smiles.substring(j,k)+tmpsmiles.substring(2,tmpsmiles.length());
                                i++;
                                continue;
                            }
                            else{
                                tmpsmiles = tmpsmiles.substring(0,1)+smiles.substring(k,i+1) + smiles.substring(j,k)+tmpsmiles.substring(1,tmpsmiles.length());
                                i++;
                                continue;
                            }
                        }
                    }
                    if(tmpc1.equals('.')){
                        int k = i+1;
                        while(!tmpc1.equals(']')){
                            i++;
                            tmpc1 = smiles.charAt(i);
                        }
                        if(j==0){
                            tmpsmiles = smiles.substring(k,i+1) + "." + smiles.substring(j,k);
                            i++;
                            continue;
                        }
                        if(j==1){
                            tmpsmiles = smiles.substring(k,i+1) + "." + smiles.substring(j,k) + smiles.substring(0,j);
                            i++;
                            continue;
                        }
                        if(j>1){
                            Character cx = tmpsmiles.charAt(j-1);
                            if(cx>'a'&&cx<'z'||Character.isDigit(cx)){
                                tmpsmiles = tmpsmiles.substring(0,2)+smiles.substring(k,i+1) + "." + smiles.substring(j,k-1)+tmpsmiles.substring(2,tmpsmiles.length());
                                i++;
                                continue;
                            }
                            else{
                                tmpsmiles = tmpsmiles.substring(0,1)+smiles.substring(k,i+1) + "." + smiles.substring(j,k-1)+tmpsmiles.substring(1,tmpsmiles.length());
                                i++;
                                continue;
                            }
                        }
                    }
                    else{
                        if(j==0){
                            tmpsmiles = smiles.substring(0,i);
                            continue;
                        }
                        if(j==1){
                            tmpsmiles = smiles.substring(0,i);
                            continue;
                        }
                        if(j>1){
                            Character cx = tmpsmiles.charAt(j-1);
                            if(cx>'a'&&cx<'z'||Character.isDigit(cx)){
                                tmpsmiles = tmpsmiles.substring(0,2)+smiles.substring(j,i)+tmpsmiles.substring(2,tmpsmiles.length());
                                continue;
                            }
                            else{
                                tmpsmiles = tmpsmiles.substring(0,1)+smiles.substring(j,i)+tmpsmiles.substring(1,tmpsmiles.length());
                                continue;
                            }
                        }
                    }
                }
                if(tmpc1>'A'&&tmpc1<'Z'){
                    //System.out.println(tmpc2);
                    if(tmpc2>'a'&&tmpc2<'z'){

                        tmpsmiles = "" + tmpc1 + tmpc2 + tmpsmiles;
                        //System.out.println(tmpsmiles);
                        i+=2;
                        continue;
                    }
                    if(Character.isDigit(tmpc2)){
                        if(Character.isDigit(tmpc3)){
                            tmpsmiles = "" + tmpc1 + tmpc2 + tmpc3 + tmpsmiles;
                            i+=3;
                            continue;
                        }
                        else{
                            tmpsmiles = "" + tmpc1 + tmpc2 + tmpsmiles;
                            i+=2;
                            continue;
                        }
                    }

                    else {
                        tmpsmiles = "" + tmpc1 + tmpsmiles;
                        i++;
                        continue;
                    }
                }
                else {
                    tmpsmiles = "" + tmpc1 + tmpsmiles;
                    i++;
                    tmpc1 = smiles.charAt(i);
                    continue;
                }
            }
        }
        return tmpsmiles;
    }

    public boolean ChiralityTest(String branch1,String branch2,String branch3,String branch4){
        String[] branchs = {branch1,branch2,branch3,branch4};

        for(int i = 0;i<4;i++){
            for(int j = i+1;j<4;j++){
                if(branchs[i].equals(branchs[j])){
                    return false;
                }
            }
        }
        //System.out.println("a");

        return true;

    }
    public boolean ChiralityTest(String branch1,String branch2,String branch3){
        if(branch1.equals(branch2)||branch1.equals(branch3)||branch2.equals(branch3)){
            return false;
        }
        else{
            return true;
        }
    }


    public void LastStep(String outfile){
        File file = new File(outfile);
        BufferedReader reader = null;
        try {
            System.out.println("以行为单位读取文件内容，一次读一整行：");
            reader = new BufferedReader(new FileReader(file));
            String tempString = null;
            int line = 1;
            // 一次读入一行，直到读入null为文件结束
            while ((tempString = reader.readLine()) != null) {
                // 显示行号
                //System.out.println("line " + line + ": " + tempString);
                line++;
            }

            String[] files = new String[line];
            for(int i = 0;i<files.length;i++){
                files[i] = "";
            }
            reader.close();
            reader = new BufferedReader(new FileReader(file));
            int count = 0;
            while ((tempString = reader.readLine()) != null) {
                // 显示行号
                //System.out.println("line " + line + ": " + tempString);
                files[count] += tempString;
                count++;
            }
            reader.close();
            count = 0;
            boolean[] tags = new boolean[line];
            for(int i = 0;i<line;i++){
                tags[i] = false;
            }
            int[] pos = new int[1000];
            int m = 0;
            while(count<line){

                while(count<file.length()&&files[count].length()<6){
                    count++;
                    if(count>=line){break;}


                }
                if(count>=line){break;}
                if(files[count].substring(0,6).equals("Number")){
                    //System.out.println(files[count]);
                    pos[m] = count;
                    m++;
                }
                count++;
            }
            int[] positions = new int[m*2];
            for(int i = 0;i<m-1;i++){
                int startpos = 0;
                //System.out.println(pos[i]);
                if(files[pos[i]+2].equals("")){continue;}
                if(files[pos[i]+2].substring(0,4).equals("Salt")){
                    for(int j = pos[i]+3;j<pos[i+1];j++){
                        if(files[j].length()<4){continue;}
                        if(files[j].substring(0,4).equals("$$$$")){
                            startpos = j+1;
                            break;
                        }
                    }
                }
                else{
                    startpos = pos[i]+2;
                }
                positions[2*i] = startpos;
                positions[2*i+1] = pos[i+1]-2;
                if(i==m-2){
                    if(files[pos[i+1]+2].substring(0,4).equals("Salt")){
                        for(int j = pos[i]+3;j<files.length;j++){
                            if(files[j].length()<4){continue;}
                            if(files[j].substring(0,4).equals("$$$$")){
                                startpos = j+1;
                                break;
                            }
                        }
                    }
                    else{
                        startpos = pos[i+1]+2;
                    }
                    positions[2*i+2] = startpos;
                    positions[2*i+3] = files.length-2;
                }
            }
            for(int i = 0;i<m;i++){
                for(int j = positions[i*2];j<=positions[i*2+1];j++){
                    for(int k = j+1;k<=positions[i*2+1];k++){
                        if(files[j].equals(files[k])){
                            tags[k] = true;
                        }
                    }
                }
            }
            for(int i = 0;i< files.length;i++){
                System.out.println(tags[i]+""+i);

            }
            try {
                //如果文件存在，则追加内容；如果文件不存在，则创建文件
                File f=new File("D://outer4.sdf");
                FileWriter fw = new FileWriter(f, true);
                PrintWriter pw = new PrintWriter(fw);
                FileWriter fw2 = new FileWriter(f);
                PrintWriter pw2 = new PrintWriter(fw2);
                pw2.println("");
                pw2.flush();
                fw2.flush();
                pw2.close();
                fw2.close();
                count = 0;
                while(count<line){
                    if(tags[count]){
                        //System.out.println(count);
                        count++;
                        continue;
                    }
                    pw.println(files[count]);
                    pw.flush();
                    fw.flush();
                    count++;
                }
                pw.close();
                fw.close();
            } catch (IOException e) {
                e.printStackTrace();
            }




        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e1) {
                }
            }
        }



    }


    public boolean ChiralityCheck2(String smiles,String Presmiles){
        /*
        * 递归分支结果的使用问题
        * 环状结构处理没加
        *
        *
        *
        * */
        //System.out.println(smiles);
        int i = 0;
        boolean flag = false;
        boolean ThisB = false;
        while(i<smiles.length()){
            //System.out.println(i);
            //System.out.println(smiles);
            flag = false;
            boolean hasUSR = false;
            if(i<smiles.length()-1){

                Character ctmp = smiles.charAt(i);
                Character ctmp2 = smiles.charAt(i+1);
                //if(Character.isDigit(ctmp)){}
                if(!ctmp.equals('C')){
                    //System.out.println("haha:"+i);
                    if(Character.isDigit(ctmp)&&i>1&&i<smiles.length()-2){
                        Character tmpx = smiles.charAt(i-2);
                        Character tmpy = smiles.charAt(i+1);
                        Character tmpz = smiles.charAt(i-1);
                        if(!tmpx.equals('=')&&!tmpy.equals('=')&&tmpz.equals('C')){
                            flag = true;
                            ThisB = true;
                            pw.println("Chirality enumeration: "+i + " in " + smiles);
                            pw.flush();
                            try {
                                fw.flush();
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        }
                    }
                    if(ctmp.equals('(')&&i!=0){
                        Character cx = smiles.charAt(i-1);
                        if(cx.equals('N')){
                            i++;continue;
                        }

                    }
                    }
                if(ctmp<='z'&&ctmp>='a'||ctmp2<='z'&&ctmp2>='a'){i++;continue;}
            }
            if(i==smiles.length()-1){
                Character ctmp = smiles.charAt(i);
                if(!ctmp.equals('C')&&!Character.isDigit(ctmp)){i++;continue;}
                if(ctmp<='z'&&ctmp>='a'){i++;continue;}
            }
            //System.out.println("hoho:"+i+smiles);
            int[] bpos = SubBranchFind(smiles,i);
            if(bpos[0]==-1){
                i++;
                continue;

            }
            //System.out.println("heihei:"+bpos[bpos.length-1]);
            if(bpos.length==1){
                i++;
                continue;
            }
            //System.out.println(bpos.length+" "+i);
            if(bpos.length==2){
                //boolean flag = false;
                String[] branchs = new String[3];
                String[] Prebranchs = new String[2];
                branchs[0] = ReservePreBranch(smiles,bpos[0]);
                //System.out.println("A");
                branchs[0] = branchs[0] + Presmiles;
                if(branchs[0].equals("")){
                    flag = false;
                }
                //branchs[0] = RingStandard(branchs[0]);
                branchs[1] = CatchBranch(smiles,bpos[0]);
                //branchs[1] = RingStandard(branchs[1]);
                branchs[2] = CatchBranch(smiles,bpos[1]);
                //branchs[2] = RingStandard(branchs[2]);
                /*
                if(branchs[2].equals("")||branchs[0].equals("")||branchs[1].equals("")){
                    i++;
                    continue;
                }
                */
                Character tmp = smiles.charAt(bpos[0]-1);
                if(Character.isDigit(tmp)){
                    Prebranchs[0] =  "C" + tmp + "(" + branchs[0] + ")" + branchs[2];
                    Prebranchs[1] =  "C" + tmp + "(" + branchs[0] + ")" + branchs[1];
                }
                else{
                    Prebranchs[0] =  tmp + "(" + branchs[0] + ")" + branchs[2];
                    Prebranchs[1] =  tmp + "(" + branchs[0] + ")" + branchs[1];
                }

                if(branchs[0].equals("")){branchs[0] = "H";}
                if(branchs[1].equals("")){branchs[1] = "H";}
                if(branchs[2].equals("")){branchs[2] = "H";}


                if(branchs[0].substring(0,1).equals("=")){
                    flag = false;
                    ChiralityCheck2(branchs[1],Prebranchs[0]);
                    ChiralityCheck2(branchs[2],Prebranchs[1]);
                    i = bpos[1]+branchs[2].length()+1;
                    continue;
                }
                if(branchs[1].substring(0,1).equals("=")){
                    flag = false;
                    ChiralityCheck2(branchs[1],Prebranchs[0]);
                    ChiralityCheck2(branchs[2],Prebranchs[1]);
                    i = bpos[1]+branchs[2].length()+1;
                    continue;
                }
                if(branchs[2].substring(0,1).equals("=")){
                    flag = false;
                    ChiralityCheck2(branchs[1],Prebranchs[0]);
                    ChiralityCheck2(branchs[2],Prebranchs[1]);
                    i = bpos[1]+branchs[2].length()+1;
                    continue;
                }

                if(!IsUSRIn){
                    if(!IsASRING(branchs[0])){flag = false;hasUSR = true;}
                    if(!IsASRING(branchs[1])){flag = false;hasUSR = true;}
                    if(!IsASRING(branchs[2])){flag = false;hasUSR = true;}
                }
                if(IsUSRIn){
                    if(!IsASRING(branchs[0])){
                        flag = true;
                        pw.println("Chirality enumeration: "+i + " in " + smiles);
                        pw.flush();
                        try {
                            fw.flush();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        //ChiralityCheck2(branchs[0]);
                    }
                    else if(!IsASRING(branchs[1])){
                        flag = true;
                        pw.println("Chirality enumeration: "+i + " in " + smiles);
                        pw.flush();
                        try {
                            fw.flush();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        ChiralityCheck2(branchs[1],Prebranchs[0]);
                    }
                    else if(!IsASRING(branchs[2])){
                        flag = true;
                        pw.println("Chirality enumeration: "+i + " in " + smiles);
                        pw.flush();
                        try {
                            fw.flush();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        ChiralityCheck2(branchs[2],Prebranchs[1]);
                    }
                }
                int i1 = branchs[0].length();
                int i2 = branchs[1].length();
                int i3 = branchs[2].length();
                /*
                if(ChiralityCheck2(branchs[0])){
                    flag = true;
                    //System.out.println(i + " in " + smiles);
                }
                */

                if(ChiralityCheck2(branchs[1],Prebranchs[0])){
                    flag = true;
                    pw.println("Chirality enumeration: "+i + " in " + smiles);
                    pw.flush();
                    try {
                        fw.flush();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                //System.out.println(branchs[2]);
                if(ChiralityCheck2(branchs[2],Prebranchs[1])) {
                    flag = true;
                    pw.println("Chirality enumeration: "+i + " in " + smiles);
                    pw.flush();
                    try {
                        fw.flush();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }

                if(!flag&&!hasUSR){
                    if(ChiralityTest(RingStandard(branchs[0]), RingStandard(branchs[1]), RingStandard(branchs[2]))){
                        //System.out.println(i + " in " + smiles);
                        pw.println("Chirality enumeration: "+i + " in " + smiles);
                        pw.flush();
                        try {
                            fw.flush();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        flag = true;
                        //return true;
                    }
                }
                if(flag){
                    ThisB = true;
                }
                i = bpos[1]+branchs[2].length()+1;
            }
            if(bpos.length==3) {
                //boolean flag = false;
                String[] branchs = new String[4];
                String[] Prebranchs = new String[3];
                branchs[0] = ReservePreBranch(smiles, bpos[0]);
                branchs[0] = branchs[0] + Presmiles;
                if(branchs[0].equals("")){flag = false;}
                //branchs[0] = RingStandard(branchs[0]);
                branchs[1] = CatchBranch(smiles, bpos[0]);
                //branchs[1] = RingStandard(branchs[1]);
                branchs[2] = CatchBranch(smiles, bpos[1]);
                //branchs[2] = RingStandard(branchs[2]);
                branchs[3] = CatchBranch(smiles, bpos[2]);
                //branchs[3] = RingStandard(branchs[3]);
                /*
                if (ChiralityCheck2(branchs[0])) {
                    flag = true;
                    //System.out.println(i + " in " + smiles);
                }

                 */
                Character tmp = smiles.charAt(bpos[0]-1);
                if(Character.isDigit(tmp)){
                    Prebranchs[0] =  "C" + tmp + "(" + branchs[0] + ")" + "(" + branchs[2] + ")"+ branchs[3];
                    Prebranchs[1] =  "C" + tmp + "(" + branchs[0] + ")" + "(" + branchs[1] + ")"+ branchs[3];
                    Prebranchs[2] =  "C" + tmp + "(" + branchs[0] + ")" + "(" + branchs[1] + ")"+ branchs[2];
                }
                else{
                    Prebranchs[0] =  tmp + "(" + branchs[0] + ")" + "(" + branchs[2] + ")" + branchs[3];
                    Prebranchs[1] =  tmp + "(" + branchs[0] + ")" + "(" + branchs[1] + ")" + branchs[3];
                    Prebranchs[2] =  tmp + "(" + branchs[0] + ")" + "(" + branchs[1] + ")" + branchs[2];
                }
                if(branchs[3].equals("")||branchs[2].equals("")||branchs[0].equals("")||branchs[1].equals("")){
                    i++;
                    continue;
                }
                if(branchs[0].substring(0,1).equals("=")){
                    flag = false;
                    ChiralityCheck2(branchs[1],Prebranchs[0]);
                    ChiralityCheck2(branchs[2],Prebranchs[1]);
                    ChiralityCheck2(branchs[3],Prebranchs[2]);
                    i = bpos[2]+branchs[3].length()+1;
                    continue;
                }
                if(branchs[1].substring(0,1).equals("=")){
                    flag = false;
                    ChiralityCheck2(branchs[1],Prebranchs[0]);
                    ChiralityCheck2(branchs[2],Prebranchs[1]);
                    ChiralityCheck2(branchs[3],Prebranchs[2]);
                    i = bpos[2]+branchs[3].length()+1;
                    continue;
                }
                if(branchs[2].substring(0,1).equals("=")){
                    flag = false;
                    ChiralityCheck2(branchs[1],Prebranchs[0]);
                    ChiralityCheck2(branchs[2],Prebranchs[1]);
                    ChiralityCheck2(branchs[3],Prebranchs[2]);
                    i = bpos[2]+branchs[3].length()+1;
                    continue;
                }
                if(branchs[3].substring(0,1).equals("=")){
                    flag = false;
                    ChiralityCheck2(branchs[1],Prebranchs[0]);
                    ChiralityCheck2(branchs[2],Prebranchs[1]);
                    ChiralityCheck2(branchs[3],Prebranchs[2]);
                    i = bpos[2]+branchs[3].length()+1;
                    continue;
                }

                if(!IsUSRIn){
                    if(!IsASRING(branchs[0])){flag = false;hasUSR = true;}
                    if(!IsASRING(branchs[1])){flag = false;hasUSR = true;}
                    if(!IsASRING(branchs[2])){flag = false;hasUSR = true;}
                    if(!IsASRING(branchs[3])){flag = false;hasUSR = true;}
                }
                if(IsUSRIn){
                    if(!IsASRING(branchs[0])){
                        flag = true;
                        pw.println("Chirality enumeration: "+i + " in " + smiles);
                        pw.flush();
                        try {
                            fw.flush();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        //ChiralityCheck2(branchs[0]);
                    }
                    else if(!IsASRING(branchs[1])){
                        flag = true;
                        pw.println("Chirality enumeration: "+i + " in " + smiles);
                        pw.flush();
                        try {
                            fw.flush();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        ChiralityCheck2(branchs[1],Prebranchs[0]);
                    }
                    else if(!IsASRING(branchs[2])){
                        flag = true;
                        pw.println("Chirality enumeration: "+i + " in " + smiles);
                        pw.flush();
                        try {
                            fw.flush();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        ChiralityCheck2(branchs[2],Prebranchs[1]);
                    }
                    else if(!IsASRING(branchs[3])){
                        flag = true;
                        pw.println("Chirality enumeration: "+i + " in " + smiles);
                        pw.flush();
                        try {
                            fw.flush();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        ChiralityCheck2(branchs[3],Prebranchs[2]);
                    }
                }
                if (ChiralityCheck2(branchs[1],Prebranchs[0])) {
                    flag = true;
                    pw.println("Chirality enumeration: "+i + " in " + smiles);
                    pw.flush();
                    try {
                        fw.flush();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                if (ChiralityCheck2(branchs[2],Prebranchs[1])) {
                    flag = true;
                    pw.println("Chirality enumeration: "+i + " in " + smiles);
                    pw.flush();
                    try {
                        fw.flush();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                //System.out.println(branchs[3]);
                if (ChiralityCheck2(branchs[3],Prebranchs[2])) {
                    flag = true;
                    pw.println("Chirality enumeration: "+i + " in " + smiles);
                    pw.flush();
                    try {
                        fw.flush();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }


                if(!flag&&!hasUSR){
                    if (ChiralityTest(RingStandard(branchs[0]), RingStandard(branchs[1]), RingStandard(branchs[2]),RingStandard(branchs[3]))) {
                        //System.out.println(i + " in " + smiles);
                        pw.println("Chirality enumeration: "+i + " in " + smiles);
                        pw.flush();
                        try {
                            fw.flush();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        flag = true;
                    }
                }

                if(flag){
                    ThisB = true;
                }
                i = bpos[2]+branchs[3].length()+1;
            }
            else{
                //System.out.println("something may wrong");
                return false;
            }
        }
        return ThisB;
    }

    public boolean DounleInBranchTest(String smiles){
        for(int i = 0;i<smiles.length();i++){
            Character tmpc = smiles.charAt(i);
            if(tmpc.equals('=')){
                return true;
            }
        }
        return false;
    }

    public int[] SubDoubleBranchFind(String smiles,int startpos){
        //-1表示没有顺反，-2表示存在顺反且不需要比对
        int[] posi = new int[4];
        if(startpos<=1){
            posi[0] = 0;
            posi[1] = 0;
        }
        Character prec = ' ';
        if(startpos>0){
            prec = smiles.charAt(startpos-1);
        }

        Character nextc = smiles.charAt(startpos+1);
        if(!nextc.equals('C')&&!nextc.equals('c')){
            posi[2] = -1;
            posi[3] = -1;
        }
        else{
            int tmppos[] = SubBranchFind(smiles,startpos+2);
            if(startpos+2>=smiles.length()){
                posi[2] = -1;
                posi[3] = -1;
            }
            else if(tmppos.length==1){
                Character tmpc = smiles.charAt(startpos+2);
                if(tmpc.equals(')')){
                    posi[2] = -1;
                    posi[3] = -1;
                }
                else{
                    posi[2] = startpos+1;
                    posi[3] = -2;
                }
            }
            if(tmppos.length==2){
                posi[2] = tmppos[0];
                posi[3] = tmppos[1];
            }
        }
        if(!prec.equals(')')){
            if(startpos<=1){
                posi[0] = 0;
                posi[1] = 0;
            }
            else{
                posi[0] = 0;
                posi[1] = -2;

            }
        }
        else{
            int balance = -1;
            int i = startpos-2;
            Character tmpc = smiles.charAt(i);
            while(balance!=0){
                if(tmpc.equals('('))
                {
                    balance++;
                }
                if(tmpc.equals(')'))
                {
                    balance--;
                }
                i--;
                tmpc = smiles.charAt(i);
            }
            if(i==0){
                posi[0] = -2;
                posi[1] = 1;
                System.out.println("When testing Cis,smiles has some problems.");
            }
            else{
                posi[0] = 0;
                posi[1] = i+1;
            }
        }
        return posi;
    }

    public int[] FindmyBranch(String smiles,int startpos){
        int[] mypos = new int[2];
        mypos[0] = 0;
        mypos[1] = smiles.length()-1;
        if(startpos>=smiles.length()/2){
            int i = startpos+1;
            Character tmpi = smiles.charAt(i);
            int balance = 0;
            while(balance>=0&&i<smiles.length()){
                tmpi = smiles.charAt(i);
                if(tmpi.equals('(')){
                    balance ++;
                }
                if(tmpi.equals(')')){
                    balance --;
                }
                i++;

            }
            if(balance<0){
                mypos[1] = i-2;
                i = startpos-1;
                tmpi = smiles.charAt(i);
                balance = 0;
                while(balance<=0&&i>=0){
                    tmpi = smiles.charAt(i);
                    if(tmpi.equals('(')){
                        balance ++;
                    }
                    if(tmpi.equals(')')){
                        balance --;
                    }
                    i--;

                }
                if(balance>0){
                    mypos[0] = i+2;
                }
            }
        }
        else {
            int i = startpos-1;
            Character tmpi = smiles.charAt(i);
            int balance = 0;
            while(balance<=0&&i>=0){
                tmpi = smiles.charAt(i);
                if(tmpi.equals('(')){
                    balance ++;
                }
                if(tmpi.equals(')')){
                    balance --;
                }
                i--;

            }
            if(balance>0){
                mypos[0] = i+2;
                i = startpos+1;
                tmpi = smiles.charAt(i);
                balance = 0;
                while(balance>=0&&i<smiles.length()){
                    tmpi = smiles.charAt(i);
                    if(tmpi.equals('(')){
                        balance ++;
                    }
                    if(tmpi.equals(')')){
                        balance --;
                    }
                    i++;

                }
                if(balance<0){
                    mypos[1] = i-2;
                }
            }
        }
        return mypos;
    }

    public String CutBranch(String smiles,int startpos,int endpos){
        String tmp = "";
        for(int i = startpos;i<=endpos;i++){
            Character ctmp = smiles.charAt(i);
            tmp += ctmp;
        }
        return tmp;
    }

    public int[] FindRingsRange(String smiles){
        Character cx;
        int count = 0;
        int[] range = new int[200];
        for (int i = 0;i<200;i++){
            range[i] = -1;
        }
        for(int i = 0;i<smiles.length()-1;i++){
            cx = smiles.charAt(i);
            if(Character.isDigit(cx)){
                boolean flag = false;
                for(int n = 0;n<range.length;n++){
                    if(range[n]==i){flag = true;break;}
                }
                if(flag){continue;}
                range[count] = i;
                count++;
                Character cy;
                boolean doub = false;
                for(int j = i+1;j<smiles.length();j++){
                    cy = smiles.charAt(j);
                    if(cx.equals(cy)){
                        range[count] = j;
                        count++;
                        doub = true;
                        //System.out.println(i+""+j);
                    }
                    else if(!doub&&j==smiles.length()-1){
                        //System.out.println(i+""+j);
                        count--;
                        range[count] = -1;
                    }
                }
            }
        }
        int num = 0;
        for(int i = 0;i<200;i++){
            if(range[i]==-1){
                break;
            }
            num++;
        }
        int[] sam = new int[num];
        for(int i = 0;i<num;i++){
            sam[i] = range[i];
        }
        return sam;
    }

    public int[] FindBenRing(String smiles){
        int[] range = new int[200];
        for(int i = 0;i<200;i++){
            range[i] = -1;
        }
        int count = 1;
        int rcount = 0;
        int i = 0;
        boolean flag = true;
        while(flag){
            //System.out.println(count);
            i = 1;
            Character cx = smiles.charAt(i);
            while(i<smiles.length()){

                cx = smiles.charAt(i);
                //System.out.println(cx+"A");
                if((cx+"").equals(count+"")){
                    //System.out.println(i);
                    int j = i+1;
                    Character cy = smiles.charAt(j);

                    while(!(cy+"").equals(count+"")&&j<smiles.length()){
                        //System.out.println(i+"B");
                        j++;
                        cy = smiles.charAt(j);
                    }
                    if(cx.equals(cy)){
                        //System.out.println("C++");
                        cx = smiles.charAt(i-1);
                        cy = smiles.charAt(j-1);
                        if(!cx.equals('C')||!cy.equals('C')){
                            count++;
                            break;
                        }
                        int eqcou = 0;
                        int m = i;
                        int Ccount = 2;
                        while(m<j){
                            //System.out.println(m+j);
                            cx = smiles.charAt(m);
                            //System.out.println(cx+" "+i+" "+j+" "+Ccount + " " + eqcou);
                            if(cx.equals('(')){
                                int bal = -1;
                                m++;
                                while(bal<0){
                                    cx = smiles.charAt(m);
                                    //System.out.println(cx+" "+i+" "+j);
                                    if(cx.equals('(')){
                                        bal--;
                                    }
                                    if(cx.equals(')')){
                                        bal++;
                                    }
                                    m++;
                                }
                                cx = smiles.charAt(m);
                                continue;
                            }
                            if(cx.equals('=')){
                                eqcou++;
                                cy = smiles.charAt(m+1);
                                if(!cy.equals('C')){
                                    eqcou+=100;
                                }
                                if(Ccount!=2){
                                    eqcou+=100;
                                }
                                if(Ccount==2){
                                    Ccount = 0;
                                }

                            }
                            if(cx.equals('C')){
                                Ccount++;
                            }

                            m++;
                        }
                        //System.out.println(eqcou+"D");
                        if(eqcou==3){
                            range[rcount*2] = i;
                            range[rcount*2+1] = j;
                            //System.out.println(range[rcount*2+1]);
                            rcount++;
                            count++;
                            break;
                        }
                        else{
                            count++;
                            break;
                        }
                    }
                    else {
                        System.out.println("error!");
                        count++;
                        break;
                    }
                }
                i++;
                if(i==smiles.length()-1){
                    flag = false;
                }

            }
        }
        int num = 0;
        for(int m = 0;m<200;m++){
            if(range[m]==-1){
                count++;
                break;
            }
            num++;
        }
        int[] sam = new int[num];
        for(int m = 0;m<num;m++){
            sam[m] = range[m];
        }
        return sam;

    }

    public boolean DoubleNextCarbonTest(String smiles,int startpos){
        Character cx = smiles.charAt(startpos-1);
        if(Character.isDigit(cx)){
            cx = smiles.charAt(startpos-2);
            if(Character.isDigit(cx)){
                cx = smiles.charAt(startpos-3);
                if(!cx.equals('C')){
                    return false;
                }
            }
            else if(!cx.equals('C')){
                return false;
            }
        }
        else if(!cx.equals('C')){
            return false;
        }
        cx = smiles.charAt(startpos+1);
        if(Character.isDigit(cx)){
            cx = smiles.charAt(startpos+2);
            if(Character.isDigit(cx)){
                cx = smiles.charAt(startpos+3);
                if(!cx.equals('C')){
                    return false;
                }
            }
            else if(!cx.equals('C')){
                return false;
            }
        }
        else if(!cx.equals('C')){
            return false;
        }
        return true;
    }

    public boolean DoubleNextCarbonTest2(String smiles,int startpos){
        Character cx = smiles.charAt(startpos-1);
        if(Character.isDigit(cx)){
            return false;
        }
        cx = smiles.charAt(startpos+1);
        if(Character.isDigit(cx)){
            return false;
        }
        if(startpos+2<smiles.length()){
            cx = smiles.charAt(startpos+2);
            if(Character.isDigit(cx)){
                return false;
            }

        }
        return true;
    }

    public boolean CheckDoubleInRing(int pos,int[] range){
        for(int i = 0;i<range.length/2;i++){
            if(pos>=range[i*2]&&pos<=range[i*2+1]){
                return true;
            }
        }
        return false;
    }

    public boolean CisTest(String smiles,int startpos){

        if(!DoubleNextCarbonTest(smiles,startpos)){
            return false;
        }
        if(!DoubleNextCarbonTest2(smiles,startpos)){
            return false;
        }
        int[] range = FindBenRing(smiles);
        if(CheckDoubleInRing(startpos,range)){
            return false;
        }


        int[] tmppos = FindmyBranch(smiles,startpos);
        String tmpsmile = CutBranch(smiles,tmppos[0],tmppos[1]);
        int i = startpos - tmppos[0];
        int[] sdb = SubDoubleBranchFind(tmpsmile,i);
        boolean flagleft = false;
        boolean flagright = false;
        if(sdb[0]==-1&&sdb[1]==-1){
            flagleft = false;
            return false;
        }
        if(sdb[2]==-1&&sdb[3]==-1){
            flagright = false;
            return false;
        }
        if(sdb[0]==-2||sdb[1]==-2){
            flagleft = true;
        }
        if(sdb[2]==-2||sdb[3]==-2){
            flagright = true;
        }
        if(sdb[0]>=0&&sdb[1]>=0&&tmppos[0]>1){
            String branch1 = ReservePreBranch(tmpsmile,sdb[1]) + smiles.substring(tmppos[0]-2,tmppos[0]-1) + "(" + ReservePreBranch(smiles,tmppos[0]-1) + ")" + smiles.substring(tmppos[1]+2,smiles.length());
            String branch2 = CatchBranch(tmpsmile,sdb[1]);
            if(!branch1.equals(branch2)){
                flagleft = true;
            }
        }
        if(sdb[2]>=0&&sdb[3]>=0){
            String branch1 = CatchBranch(tmpsmile,sdb[2]);
            String branch2 = CatchBranch(tmpsmile,sdb[3]);
            if(!branch1.equals(branch2)){
                flagright = true;
            }
        }
        if(flagleft&&flagright){
            return true;
        }
        return false;
    }

    public String CutHeadSalt(String smiles){
        System.out.println(smiles);
        Character c0 = smiles.charAt(0);
        Character c1 = smiles.charAt(1);
        Character c2 = smiles.charAt(2);
        if(c1.equals('.')){
            return smiles.substring(1,smiles.length());
        }
        if(c2.equals('.')){
            return smiles.substring(2,smiles.length());
        }
        if(!c0.equals('[')){
            System.out.println("No [] on its head!");
            return smiles;
        }
        else {
            String tmpsmiles = "";
            boolean start = false;
            for (int i = 1; i < smiles.length(); i++) {
                Character ct = smiles.charAt(i);
                if (start) {
                    tmpsmiles += ct;
                } else {
                    if (ct.equals(']')&&!start) {
                        Character c = smiles.charAt(i+1);
                        if(!c.equals('[')){
                            i = 0;
                            ct = smiles.charAt(0);
                            tmpsmiles += ct;
                        }
                        start = true;
                    }
                }
            }
            return tmpsmiles;
        }
    }

    public String CatchSaltBranch(String smiles,int startpos){
        //startpos从[开始
        int leftcount = 1;
        int rightcount = 0;
        int i = startpos+1;
        String tmpsmiles = "";
        if(startpos==smiles.length()-1){
            return "" + smiles.charAt(i-1);
        }

        Character tmps = smiles.charAt(i);
        Character tmps2 = smiles.charAt(startpos);
        if(!tmps2.equals('[')&&!tmps2.equals(']')){
            tmpsmiles += tmps2;
        }
        while(leftcount!=rightcount&&smiles.length()>i){
            tmps = smiles.charAt(i);
            if(tmps.equals('[')){
                leftcount++;
            }
            if(tmps.equals(']')){
                rightcount++;
            }
            if(leftcount!=rightcount) {
                tmpsmiles += tmps;
            }
            i++;
        }
        return  tmpsmiles;
    }

    public int CatchEndPos(String smiles,int startpos){
        //startpos从[开始
        int leftcount = 1;
        int rightcount = 0;
        int i = startpos+1;
        if(startpos==smiles.length()-1){
            return startpos-1;
        }

        Character tmps = smiles.charAt(i);
        Character tmps2 = smiles.charAt(startpos);
        if(!tmps2.equals('[')&&!tmps2.equals(']')){

        }
        while(leftcount!=rightcount&&smiles.length()>i){
            tmps = smiles.charAt(i);
            if(tmps.equals('[')){
                leftcount++;
            }
            if(tmps.equals(']')){
                rightcount++;
            }
            if(leftcount!=rightcount) {

            }
            i++;
        }
        return  --i;
    }

    public String RemoveSalt(String smiles){
        String smilesch = CutHeadSalt(smiles);
        int i = 1;
        String tmpsmiles = "";
        Character tmpc = smilesch.charAt(0);
        Character ct = smilesch.charAt(0);
        if(tmpc.equals('.')){
            i = 2;
            ct = smilesch.charAt(1);
        }
        tmpsmiles += ct;
        ct = smilesch.charAt(i);
        while(i<smilesch.length()){
            if(!ct.equals('[')){
                tmpsmiles += ct;
                i++;
                if(i<smilesch.length()){
                    ct = smilesch.charAt(i);
                }

            }
            else if(ct.equals('[')){
                //System.out.println(i);
                int j = CatchEndPos(smilesch,i);
                Character cx = smilesch.charAt(i-1);
                if(j==smilesch.length()-1){
                    tmpc = smilesch.charAt(i-1);
                    if(tmpc.equals('.')){
                        tmpsmiles = tmpsmiles.substring(0,tmpsmiles.length()-1);
                        return tmpsmiles;
                    }
                    else {return tmpsmiles;}
                }
                if(cx.equals('(')){
                    cx = smilesch.charAt(j+1);
                    if(cx.equals('[')){
                        tmpsmiles += '[';
                        tmpsmiles += CatchSaltBranch(smilesch,i);
                        tmpsmiles += ']';
                        j = CatchEndPos(smilesch,j+1);
                        i = j+1;
                        ct = smilesch.charAt(i);
                    }
                    else{
                        tmpsmiles += ct;
                        i++;
                        ct = smilesch.charAt(i);
                    }
                }
                else{

                    cx = smilesch.charAt(j+1);
                    if(cx.equals('[')){
                        tmpsmiles += '[';
                        tmpsmiles += CatchSaltBranch(smilesch,i);
                        tmpsmiles += ']';
                        if(tmpsmiles.length()>=smilesch.length()-j-1){
                            return tmpsmiles;
                        }
                        else {
                            i = j+2;
                            tmpsmiles = ""+'[';
                            ct = smilesch.charAt(i);
                        }
                    }
                    else{
                        i++;
                        tmpsmiles += ct;
                        ct = smilesch.charAt(i);
                    }


                }
            }
        }
        return tmpsmiles;
    }


    public void AllCheck(boolean ifout) throws IOException, CDKException {
        //System.out.println(ifout);
        if(smiles.equals("")){return;}
        if(ifout){return;}
        pw.println(smiles);
        pw.flush();
        try {
            fw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

        //System.out.println(smiles);
        String mysmiles = RemoveSalt(smiles);
        if(!mysmiles.equals(smiles)){
            pw.println("Salt processing:"+smiles+" has removed little structures");
            //pw.println();
            /*
            try {
                SmilesParser sp  = new SmilesParser(SilentChemObjectBuilder.getInstance());
                IAtomContainer m   = sp.parseSmiles(mysmiles);
                //IAtomContainerSet cs = new AtomContainerSet();
                //cs.addAtomContainer(m);
                BufferedWriter bw = new BufferedWriter(fw);
                SDFWriter sdfw = new SDFWriter(bw);
                sdfw.write(m);
                bw.flush();
                //pw.println(mysmiles);
                //pw.flush();
            } catch (InvalidSmilesException e) {
                System.err.println(e.getMessage());
            }
            */
        }
        //System.out.println("1");
        //System.out.println(mysmiles);
        if(!ifout){
            InPutSdf(mysmiles);
        }

        ChiralityCheck2(mysmiles,"");

        //System.out.println("2");

        for (int i = 0; i<mysmiles.length();i++){

            Character c = mysmiles.charAt(i);
            if(c.equals('=')){
                //System.out.println(i);
                if(CisTest(mysmiles,i)){

                    pw.println("Cis-Trans enumeration: "+i+" in "+mysmiles);
                    pw.flush();
                    try {
                        fw.flush();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }

            }
        }
        BondCountTest();
        SpecialSelect();
        TautomerTest();

    }
    public void GoNext() throws IOException {
        pw.close();
        fw.close();
        pw2.close();
        fw2.close();
    }


}



