package edu.illinois.cs.cs125.spring2019.mp4;

import android.annotation.SuppressLint;
import android.app.Activity;
import android.content.SharedPreferences;
import android.os.Bundle;
import android.util.Log;
import android.view.View;
import android.webkit.WebView;
import android.webkit.WebViewClient;
import android.widget.Switch;
import android.widget.TextView;

import com.google.gson.Gson;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;

import edu.illinois.cs.cs125.spring2019.mp4.lib.ChemicalElement;
import edu.illinois.cs.cs125.spring2019.mp4.lib.CyclicOrganicMoleculeBuilder;
import edu.illinois.cs.cs125.spring2019.mp4.lib.LinearOrganicMoleculeBuilder;
import edu.illinois.cs.cs125.spring2019.mp4.lib.OrganicMoleculeBuilder;

/**
 * random javadoc comment.
 */
public class MoleculeActivity extends Activity {

    /**
     * Tag for logging.
     */
    private static final String TAG = "MP4:MainActivity";

    /**
     * The name of the SharedPreferences file for the app.
     */
    private static final String PREFS_NAME = "CS125_MP4_Spring2019";

    /**
     * The name of the SharedPreferences entry that stores the last-viewed molecule index.
     */
    private static final String PREF_LAST_MOLECULE = "LastMoleculeIndex";

    /**
     * The name of the SharedPreferences entry that stores whether the user wants colored drawing.
     */
    private static final String PREF_COLOR = "UseColor";

    /**
     * The Android WebView object that displays the molecule drawn on top of a web canvas.
     */
    private WebView webView;

    /**
     * A list of MoleculeData objects where each element is a molecule that can be rendered.
     */
    private ArrayList<MoleculeData> loadedMolecules;

    /**
     * The index of the molecule in the loadedMolecules list that is currently being displayed.
     */
    private int currentMoleculeIndex = 0;

    /**
     * Flag to represent preference of whether to use CPK coloring for molecules
     * as opposed to black and white drawing.
     */
    private boolean useColor = true;

    /**
     * A SharedPreferences instance that allows access to the persistent settings for the app.
     */
    private SharedPreferences sharedPreferences;

    /**
     * Called when the activity is created. Sets up the WebView for displaying molecules.
     * <p>
     * You need to expand this function to restore the previously-viewed molecule index
     * from SharedPreferences so that the app remembers the viewed molecule across restarts.
     * @param savedInstanceState unused
     * @see <a href="https://developer.android.com/training/data-storage/shared-preferences#java">SharedPreferences</a>
     */
    @SuppressLint("SetJavascriptEnabled")
    @Override
    protected final void onCreate(final Bundle savedInstanceState) {
        // Set our layout
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_molecule);

        // Load all user-defined molecules
        loadedMolecules = loadMolecules();

        // Restore the color preference from SharedPreferences
        sharedPreferences = getSharedPreferences(PREFS_NAME, MODE_PRIVATE);
        useColor = sharedPreferences.getBoolean(PREF_COLOR, true);
        ((Switch) findViewById(R.id.useColorSwitch)).setChecked(useColor);

        // Restore the last-viewed molecule index from SharedPreferences
        currentMoleculeIndex = sharedPreferences.getInt(PREF_LAST_MOLECULE, 0);
        // Set up our WebView to display the current molecule
        webView = findViewById(R.id.webView);
        webView.getSettings().setJavaScriptEnabled(true);
        webView.getSettings().setBuiltInZoomControls(true);
        webView.setWebViewClient(new WebViewClient() {
            public void onPageFinished(final WebView unusedWebView, final String unusedString) {
                loadColorPreference(useColor);
                loadMolecule(currentMoleculeIndex);
            }
        });
        webView.loadUrl("file:///android_asset/index.html");
    }

    /**
     * Return an ArrayList of molecules that can be viewed and analyzed.
     * <p>
     * You should feel free to add your own molecules here to aid in debugging your library.
     *
     * @return a list of molecules
     */
    @SuppressWarnings("checkstyle:magicnumber")
    private ArrayList<MoleculeData> loadMolecules() {
        ArrayList<MoleculeData> newMolecules = new ArrayList<>();

        newMolecules.add(new MoleculeData(new LinearOrganicMoleculeBuilder(3)));
        newMolecules.add(new MoleculeData(new CyclicOrganicMoleculeBuilder(6)));
        newMolecules.add(new MoleculeData(new CyclicOrganicMoleculeBuilder(4)
                .addSubstituent(0, OrganicMoleculeBuilder.Substituent.createAlcohol())));
        newMolecules.add(new MoleculeData(new CyclicOrganicMoleculeBuilder(5)
                .addSubstituent(1, OrganicMoleculeBuilder.Substituent.createAlkyl(2))));
        newMolecules.add(new MoleculeData(new LinearOrganicMoleculeBuilder(4)
                .addSubstituent(3, OrganicMoleculeBuilder.Substituent.createCarbonyl())));
        newMolecules.add(new MoleculeData(new LinearOrganicMoleculeBuilder(3)
                .addSubstituent(1, OrganicMoleculeBuilder.Substituent.createAlkyl(1))));
        newMolecules.add(new MoleculeData(new LinearOrganicMoleculeBuilder(6)
                .addSubstituent(4, OrganicMoleculeBuilder.Substituent.createAlkyl(2))));
        newMolecules.add(new MoleculeData(new CyclicOrganicMoleculeBuilder(4)
                .addSubstituent(0, OrganicMoleculeBuilder.Substituent.createAlkyl(1))
                .addSubstituent(1, OrganicMoleculeBuilder.Substituent.createAlkyl(1))
                .addSubstituent(2, OrganicMoleculeBuilder.Substituent.createAlkyl(1))
                .addSubstituent(3, OrganicMoleculeBuilder.Substituent.createAlkyl(1))));
        newMolecules.add(new MoleculeData(new CyclicOrganicMoleculeBuilder(5)
                .addSubstituent(0, OrganicMoleculeBuilder.Substituent.createCarbonyl())
                .addSubstituent(1, OrganicMoleculeBuilder.Substituent.createCarbonyl())
                .addSubstituent(2, OrganicMoleculeBuilder.Substituent.createCarbonyl())
                .addSubstituent(3, OrganicMoleculeBuilder.Substituent.createCarbonyl())
                .addSubstituent(4, OrganicMoleculeBuilder.Substituent.createCarbonyl())));
        newMolecules.add(new MoleculeData(new LinearOrganicMoleculeBuilder(4)
                .addSubstituent(0, OrganicMoleculeBuilder.Substituent.createHalogen(ChemicalElement.CHLORINE))
                .addSubstituent(0, OrganicMoleculeBuilder.Substituent.createHalogen(ChemicalElement.CHLORINE))
                .addSubstituent(0, OrganicMoleculeBuilder.Substituent.createHalogen(ChemicalElement.CHLORINE))
                .addSubstituent(2, OrganicMoleculeBuilder.Substituent.createCarbonyl())));
        newMolecules.add(new MoleculeData(new LinearOrganicMoleculeBuilder(2)
                .addSubstituent(0, OrganicMoleculeBuilder.Substituent.createAlcohol())
                .addSubstituent(1, OrganicMoleculeBuilder.Substituent.createHalogen(ChemicalElement.CHLORINE))
                .addSubstituent(1, OrganicMoleculeBuilder.Substituent.createHalogen(ChemicalElement.BROMINE))));
        newMolecules.add(new MoleculeData(new LinearOrganicMoleculeBuilder(5)
                .addSubstituent(1, OrganicMoleculeBuilder.Substituent.createHalogen(ChemicalElement.BROMINE))
                .addSubstituent(2, OrganicMoleculeBuilder.Substituent.createAlkyl(1))
                .addSubstituent(3, OrganicMoleculeBuilder.Substituent.createAlkyl(3))
                .addSubstituent(4, OrganicMoleculeBuilder.Substituent.createCarbonyl())));
        newMolecules.add(new MoleculeData(new CyclicOrganicMoleculeBuilder(3)
                .addSubstituent(2, OrganicMoleculeBuilder.Substituent.createHalogen(ChemicalElement.CHLORINE))
                .addSubstituent(2, OrganicMoleculeBuilder.Substituent.createHalogen(ChemicalElement.BROMINE))
                .addSubstituent(1, OrganicMoleculeBuilder.Substituent.createHalogen(ChemicalElement.FLUORINE))
                .addSubstituent(1, OrganicMoleculeBuilder.Substituent.createAlcohol())
                .addSubstituent(0, OrganicMoleculeBuilder.Substituent.createHalogen(ChemicalElement.BROMINE))));

        // TODO: Add your own new molecules here!

        return newMolecules;
    }

    /**
     * Load and display the next molecule.
     *
     * @param unused used Android view object
     */
    public void onNextButtonClick(final View unused) {
        currentMoleculeIndex = (currentMoleculeIndex + 1) % loadedMolecules.size();
        loadMolecule(currentMoleculeIndex);
    }

    /**
     * Load and display the previous molecule.
     *
     * @param unused unused Android view object
     */
    public void onPrevButtonClick(final View unused) {
        currentMoleculeIndex = (currentMoleculeIndex - 1 + loadedMolecules.size()) % loadedMolecules.size();
        loadMolecule(currentMoleculeIndex);
    }

    /**
     * Toggles whether or not the molecules on screen are colored. Called when the user toggles the Switch.
     *
     * @param unused unused Android view object
     */
    public void onColorPreferenceToggle(final View unused) {
        useColor = !useColor;

        // Write the changed setting to SharedPreferences
        SharedPreferences.Editor editor = sharedPreferences.edit();
        editor.putBoolean(PREF_COLOR, useColor);
        editor.apply();

        loadColorPreference(useColor);
        loadMolecule(currentMoleculeIndex); // Need to re-render molecule when updating preference
    }

    /**
     * Notify the sketch.js molecule drawer what the new color preference is.
     *
     * @param color if true: renders molecules in color, if false: renders molecules in black and white
     */
    private void loadColorPreference(final boolean color) {
        webView.evaluateJavascript("setColorPreference(" + color + ");",
            s -> Log.i(TAG, "evaluateJavascript completed, result: " + s));
    }

    /**
     * Loads a new molecule onto the WebView.
     *
     * @param moleculeIndex the index of the molecule to load in our molecule ArrayList
     */
    private void loadMolecule(final int moleculeIndex) {
        if (moleculeIndex < 0 || moleculeIndex > loadedMolecules.size()) {
            Log.w(TAG, "Tried to load molecule at invalid index: " + moleculeIndex);
            currentMoleculeIndex = 0; // Gets the app unstuck, in case some were removed from loadMolecules
            return;
        }

        // Record the molecule index so that it can be restored if the app is restarted
        SharedPreferences.Editor editor = sharedPreferences.edit();
        editor.putInt(PREF_LAST_MOLECULE, moleculeIndex);
        editor.apply();

        // Get the molecule from the list
        MoleculeData moleculeToLoad = loadedMolecules.get(moleculeIndex);

        // Serialize to JSON for transfer to JavaScript
        String moleculeAsJSON = new Gson().toJson(moleculeToLoad);

        // We only support API level 23 and higher, but if we wanted to support pre-19 devices,
        // we would need to check the Android version and possibly fall back to loadUrl
        webView.evaluateJavascript("loadMolecule('" + moleculeAsJSON + "');",
            s -> Log.i(TAG, "evaluateJavascript completed, result: " + s));

        // Update the UI with various pieces of information about the molecule
        final TextView formulaTextView = findViewById(R.id.text_formula);
        formulaTextView.setText(moleculeToLoad.getFormula());

        final TextView nameTextView = findViewById(R.id.text_name);
        nameTextView.setText(moleculeToLoad.getName());

        final TextView massTextView = findViewById(R.id.text_mass);
        DecimalFormat decimalFormat = new DecimalFormat("#.000");
        decimalFormat.setRoundingMode(RoundingMode.CEILING);
        massTextView.setText(String.format("%s amu", decimalFormat.format(moleculeToLoad.getMass())));
    }

    /**
     * For testing purposes only. Gets the currently selected molecule index.
     * @return the index of the molecule being viewed
     */
    public int getCurrentMoleculeIndex() {
        return currentMoleculeIndex;
    }
}
