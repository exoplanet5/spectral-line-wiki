// OCR PNG page images with the macOS Vision framework.
// usage: swift ocr_vision.swift <out.txt> <img1.png> [img2.png ...]
import Foundation
import Vision
import AppKit

let args = CommandLine.arguments
guard args.count >= 3 else { print("usage: ocr_vision <out.txt> <images...>"); exit(1) }
var out = ""
for path in args[2...] {
    guard let img = NSImage(contentsOfFile: path),
          let cg = img.cgImage(forProposedRect: nil, context: nil, hints: nil) else {
        out += "\n===== \(path) (unreadable) =====\n"; continue
    }
    let req = VNRecognizeTextRequest()
    req.recognitionLevel = .accurate
    req.usesLanguageCorrection = false
    let h = VNImageRequestHandler(cgImage: cg, options: [:])
    try? h.perform([req])
    var lines: [(CGFloat, CGFloat, String)] = []
    for obs in (req.results ?? []) {
        if let c = obs.topCandidates(1).first {
            let bb = obs.boundingBox
            lines.append((1 - bb.midY, bb.minX, c.string))
        }
    }
    // group into rows by y, then sort by x
    lines.sort { $0.0 == $1.0 ? $0.1 < $1.1 : $0.0 < $1.0 }
    var rows: [[(CGFloat, CGFloat, String)]] = []
    for l in lines {
        if var last = rows.last, let y0 = last.first?.0, abs(l.0 - y0) < 0.006 {
            last.append(l); rows[rows.count - 1] = last
        } else { rows.append([l]) }
    }
    out += "\n===== \((path as NSString).lastPathComponent) =====\n"
    for r in rows {
        out += r.sorted { $0.1 < $1.1 }.map { $0.2 }.joined(separator: "  |  ") + "\n"
    }
}
try out.write(toFile: args[1], atomically: true, encoding: .utf8)
