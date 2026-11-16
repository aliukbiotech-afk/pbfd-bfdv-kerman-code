// Builds cover_letter_PBFD.docx  (DOCX_PATH=<path to docx module> node build_cover_letter.js [outfile])
const fs = require("fs");
const path = require("path");
const D = require(process.env.DOCX_PATH || "docx");
const { Document, Packer, Paragraph, TextRun, AlignmentType, LineRuleType } = D;

const PENDING = [];
function runs(text, base = {}) {
  const parts = text.split(/(\*\*[^*]+\*\*|\*[^*]+\*|\[\[[^\]]+\]\])/g).filter(s => s !== "");
  return parts.map(s => {
    if (s.startsWith("**")) return new TextRun({ ...base, text: s.slice(2, -2), bold: true });
    if (s.startsWith("*")) return new TextRun({ ...base, text: s.slice(1, -1), italics: true });
    if (s.startsWith("[[")) { PENDING.push(s.slice(2, -2)); return new TextRun({ ...base, text: "...", highlight: "yellow" }); }
    return new TextRun({ ...base, text: s });
  });
}
const P = (t, o = {}) => new Paragraph({ children: runs(t), alignment: o.align || AlignmentType.JUSTIFIED,
  spacing: { after: o.after ?? 160, line: 300, lineRule: LineRuleType.AUTO } });

const C = [];
C.push(P("[[date]]", { align: AlignmentType.LEFT }));
C.push(P("The Editor-in-Chief\n".replace("\n", "") + "[[journal name]]", { align: AlignmentType.LEFT, after: 240 }));
C.push(P("**Submission of an original research article**", { align: AlignmentType.LEFT }));
C.push(P("**Title:** Subclinical beak and feather disease virus infection in pet psittacine birds in Kerman, south-eastern Iran: prevalence, risk factors and full-genome-guided genotyping of partial Rep sequences", { align: AlignmentType.LEFT, after: 240 }));
C.push(P("Dear Editor,", { align: AlignmentType.LEFT }));
C.push(P("We submit the above manuscript for consideration as an original research article in [[journal name]]. Beak and feather disease virus (BFDV) is the most important viral pathogen of captive and wild parrots, and apparently healthy carriers are widely believed to sustain its spread through the pet trade. Yet many molecular surveys, including earlier Iranian reports, have not consistently separated clinically normal from affected birds, and most classify the virus from short diagnostic fragments whose position in the global genome-based classification is uncertain."));
C.push(P("We screened feathers from 136 clinically healthy pet psittacines in Kerman, a province not previously studied, and detected BFDV DNA in 43 birds (31.6%; 95% CI 24.4–39.8%). Using exact and penalised methods, we found no association of positivity with sex or host group. A highly significant host effect obtained when the same data were analysed as a host × four-level band-intensity table (in the underlying thesis) rested on a sparse table that mixed infection status with PCR band intensity, which we show explicitly."));
C.push(P("The main methodological contribution is the genotyping of ten partial ORF1 (Rep) sequences by phylogenetic placement (EPA-ng) on a backbone of 429 complete BFDV genomes classified according to Shah et al. (2023), with a likelihood weight ratio for each assignment and checks against constrained and unconstrained trees. All ten sequences belong to genotype GI: eight to sub-clade GI-c, one to GI-f, and one that cannot be resolved between GI-d and GI-f. The sequences have intact Rep motifs, show no evidence of recombination, and are under strong purifying selection (dN/dS = 0.14). Close matches to strains from Italy, Saudi Arabia, Poland, Bangladesh and other Iranian provinces are consistent with introduction through the live-bird trade. We report the limits of the approach openly, including unlabelled backbone genomes and the short fragment."));
C.push(P("We believe the work will interest readers of [[journal name]] because it gives, to our knowledge, the first estimate of subclinical BFDV infection in south-eastern Iran, shows how a sparse-table analysis can mislead host-effect conclusions, and offers a reproducible way to genotype the short sequences that dominate BFDV surveillance. All sequences are in GenBank (PX634653–PX634662), and the code, notebook and results are provided in a public repository ([[repository URL and DOI]]) so that every figure and table can be regenerated."));
C.push(P("**Declarations.** The manuscript is original, is not under consideration elsewhere and has not been published, in whole or in part. All authors have read and approved the submission. The authors declare no conflict of interest. Sampling was approved by [[ethics committee and approval code]], and owners gave consent [[confirm]]. Funding: [[funding sources]]."));
C.push(P("**Suggested reviewers.** [[names, affiliations and e-mail addresses of suggested reviewers; names to exclude, if any]]."));
C.push(P("Thank you for considering our work. We look forward to your response.", { after: 240 }));
C.push(P("Yours sincerely,", { align: AlignmentType.LEFT }));
C.push(P("[[corresponding author name, affiliation, e-mail, ORCID]]", { align: AlignmentType.LEFT }));
C.push(P(""));
C.push(P("**Internal checklist (delete before sending): ** each yellow '...' above stands for one of: " + PENDING.map((t, i) => `${i + 1}) ${t.replace(/^PENDING:?\s*/, "")}`).join("; ") + ".", { align: AlignmentType.LEFT }));

const doc = new Document({
  creator: "PBFD Kerman project", title: "Cover letter",
  styles: { default: { document: { run: { font: "Times New Roman", size: 24 } } } },
  sections: [{ properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } }, children: C }]
});
const out = process.argv[2] || path.join(__dirname, "..", "..", "cover_letter_PBFD.docx");
Packer.toBuffer(doc).then(b => { fs.writeFileSync(out, b); console.log("wrote", out); });
