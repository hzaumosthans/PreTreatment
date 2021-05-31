import org.openscience.cdk.ChemObject;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import java.io.*;

public class PreTreatment {

    private static String[] smiles = new String[1000];
    private static void LastStep(String outfile) {
        File file = new File(outfile);
        int EmptyCount = 0;
        int HandCount = 0;
        int CisCount = 0;
        int SaltCount = 0;
        boolean Hf = false;
        boolean Cf = false;
        boolean Sf = false;
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
            //int[] positions = new int[m*2];

            for(int i = 0;i<m;i++){
                Hf = false;
                Cf = false;
                Sf = false;
                if(pos[i]+2>=pos[i+1]){
                    EmptyCount++;
                }
                for(int j = pos[i]+1;j<=pos[i+1]-1;j++){
                    if(!Hf){
                        if(files[j].indexOf("Chirality")!=-1){
                            HandCount++;
                            Hf = true;
                        }
                    }
                    if(!Cf){
                        if(files[j].indexOf("Cis-Trans")!=-1){
                            CisCount++;
                            Cf = true;
                        }
                    }
                    if(!Sf){
                        if(files[j].indexOf("Salt")!=-1){
                            SaltCount++;
                            Sf = true;
                        }
                    }
                    for(int k = j+1;k<=pos[i+1]-1;k++){
                        if(files[j].equals(files[k])){
                            tags[k] = true;
                        }
                    }
                }
            }
            try {
                //如果文件存在，则追加内容；如果文件不存在，则创建文件
                File f=new File(outfile);
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
                pw.println("Chirality enumerations in "+HandCount+" mols");
                pw.flush();
                fw.flush();
                pw.println("Cis-Trans enumerations in "+CisCount+" mols");
                pw.flush();
                fw.flush();
                pw.println("Salt removed in "+SaltCount+" mols");
                pw.flush();
                fw.flush();
                pw.println("The number of molecules which is empty or repetition is "+EmptyCount+".");
                pw.flush();
                fw.flush();
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
    private static boolean SameTest(String thesmiles){
        int i = 0;
        while(!smiles[i].equals("")){
            if(smiles[i].equals(thesmiles)){
                return true;
            }
            i++;
        }
        smiles[i] = thesmiles;
        return false;
    }
    public static void main(String[] args) throws IOException, CDKException {
        //String infile = args[0];
        //String outfile = args[1];
        //String SdfFile = args[2];
        //Boolean IsUSRIn = args[3];
        String infile = "D://molecules1(2).sdf";
        String outfile = "D://outer3.txt";
        String SdfFile = "D://SDF.sdf";
        Boolean IsUSRIn = true;
        IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(new File(infile)),
                DefaultChemObjectBuilder.getInstance()
        );
        int count = 0;
        IChemObject set = new ChemObject();

        for(int i = 0;i<1000;i++){
            smiles[i] = "";
        }
        while (iterator.hasNext()) {
            FileWriter fw = null;
            FileWriter fw2 = null;
            try {
                //如果文件存在，则追加内容；如果文件不存在，则创建文件
                File f=new File(outfile);
                fw = new FileWriter(f, true);
            } catch (IOException e) {
                e.printStackTrace();
            }
            try {
                //如果文件存在，则追加内容；如果文件不存在，则创建文件
                File f2=new File(SdfFile);
                fw2 = new FileWriter(f2, true);
            } catch (IOException e) {
                e.printStackTrace();
            }
            PrintWriter pw = new PrintWriter(fw);
            PrintWriter pw2 = new PrintWriter(fw2);
            System.out.println(count);
            pw.println("\n"+"Number: "+count+" :  ");
            pw.flush();
            try {
                fw.flush();
            } catch (IOException e) {
                e.printStackTrace();
            }
            pw.close();
            fw.close();
            pw2.println("\n"+"Number: "+count+" :"+"\n");
            pw2.flush();
            try {
                fw2.flush();
            } catch (IOException e) {
                e.printStackTrace();
            }
            pw2.close();
            fw2.close();
            IAtomContainer mol = iterator.next();
            IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(mol);
            System.out.println(MolecularFormulaManipulator.getString(formula));
            Mymol mm = new Mymol(mol,outfile,SdfFile,IsUSRIn);
            String thesmiles = mm.returnsmiles();
            if(!thesmiles.equals("")){
                mm.AllCheck(SameTest(thesmiles));
                mm.GoNext();
            }
            count++;
        }
        FileWriter fw = null;
        try {
            //如果文件存在，则追加内容；如果文件不存在，则创建文件
            File f=new File(outfile);
            fw = new FileWriter(f, true);
        } catch (IOException e) {
            e.printStackTrace();
        }
        PrintWriter pw = new PrintWriter(fw);
        pw.println("\n"+"The total number of molecules is "+count+".");
        pw.flush();
        try {
            fw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
        pw.close();
        fw.close();
        LastStep(outfile);
    }



}
