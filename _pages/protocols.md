---
layout: archive
title: "Protocols"
permalink: /protocols/
author_profile: true
---
<script src="https://cdn.jsdelivr.net/npm/js-sha256@0.9.0/src/sha256.min.js"></script>
<script>
    function checkPassword() {
        const obfuscatedHashParts = [
            '1003f5', '9c15', '179cf021', '9817', 
            '6e9c', 'dd0333', '936c53e8', '3921ecde',
            '31003f56', 'e9c9c', '463e77c8', '3197ef'
        ];
        const getCorrectHash = () => {
            return obfuscatedHashParts[6] + obfuscatedHashParts[7] + obfuscatedHashParts[5] + 
                   '2463e77c8' + obfuscatedHashParts[8] + obfuscatedHashParts[4] + 
                   obfuscatedHashParts[9] + obfuscatedHashParts[2] + obfuscatedHashParts[11];
        };
        const correctHash = getCorrectHash();
        const userPassword = prompt("请输入访问密码：");
        if (userPassword === null || userPassword === "") {
            return;
        }
        const userHash = sha256(userPassword);
        if (userHash === correctHash) {
            window.location.href = "../protocols/Phylogenomics.html";
        } else {
            alert("密码错误，请重试！");
        }
    }
</script>
# <a onclick="checkPassword()">Phylogenomics</a>
# <a href="../protocols/transcriptome.html" target="_blank">Transcriptome</a>
# <a href="../protocols/Molecular Convergence (conv_cal).html" target="_blank">Molecular Convergence (conv_cal)</a>
# <a href="../protocols/Relaxed Selection Test.html" target="_blank">Relaxed Selection Test</a>
